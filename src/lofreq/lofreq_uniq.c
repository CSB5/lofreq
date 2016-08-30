/* -*- c-file-style: "k&r"; indent-tabs-mode: nil; -*- */
/*********************************************************************
* The MIT License (MIT)
* 
* Copyright (c) 2013,2014 Genome Institute of Singapore
* 
* Permission is hereby granted, free of charge, to any person
* obtaining a copy of this software and associated documentation files
* (the "Software"), to deal in the Software without restriction,
* including without limitation the rights to use, copy, modify, merge,
* publish, distribute, sublicense, and/or sell copies of the Software,
* and to permit persons to whom the Software is furnished to do so,
* subject to the following conditions:
* 
* The above copyright notice and this permission notice shall be
* included in all copies or substantial portions of the Software.
* 
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
* MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
* NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
* BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
* ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
* CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
* SOFTWARE.
*
************************************************************************/


/*
 * Notes about a potentially more refined analysis: Currently we take the
 * frequency for one of the SNVs as a given. Ideally we would
 * integrate over the various values possible rather then just taking the
 * maximum-likelihood value.
 *
 * Frequency estimate from first sample snv call: p^hat = k_1/n_1.
 * Current test: P_bin(n_2, p^hat) (X<=k_2).
 *
 * Better: Sum over k=0 to n_1 ( P_bin(n_2, k/n^1) (X<=k) ) * P(k).
 *
 * P(k) is prior from Binomial proportion confidence distribution
 *
 * TODO: find Binomial proportion confidence distribution to get prior
 *
 */

#include <math.h>
#include <stdio.h>
#include <unistd.h>
#include <ctype.h>
#include <string.h>
#include <errno.h>
#include <assert.h>
#include <stdarg.h>
#include <limits.h>
#include <float.h>
#include <getopt.h>
#include <stdlib.h>

/* lofreq includes */
#include "vcf.h"
#include "utils.h"
#include "log.h"
#include "binom.h"
#include "plp.h"
#include "defaults.h"
#include "snpcaller.h"
#include "multtest.h"

#if 1
#define MYNAME "lofreq uniq"
#else
#define MYNAME PACKAGE
#endif

#define DEFAULT_UNI_FREQ -1.0

#define BUF_SIZE 1<<16

#define FILTER_ID_STRSIZE 64
#define FILTER_STRSIZE 128

const char *uniq_flag = "UNIQ";

const char *uniq_phred_tag = "UQ";


typedef struct {
     int thresh;/* use if > 0; otherwise use multiple testing correction that's if >0 */
     int mtc_type;/* holm; holmbonf; fdr; none */
     double alpha;
     long int ntests;
     char id[FILTER_ID_STRSIZE];
} uniq_filter_t;


typedef struct {
     float uni_freq;
     vcf_file_t vcf_out;
     vcf_file_t vcf_in;
     int use_det_lim;
     int output_all; /* catch! doesn't actually work if there's no coverage in BAM because mpileup will skip target function */
     uniq_filter_t uniq_filter;
     /* changing per pos: the var to test */
     var_t *var;
} uniq_conf_t;





int
uniq_phred_from_var(var_t *var) {
     char *uq_char = NULL;
     if ( ! vcf_var_has_info_key(&uq_char, var, uniq_phred_tag)) {
          /* missing because no coverage or other reasons. not unique anyway */
          return 0;
     } else {
          int uq = (int) strtol(uq_char, (char **)NULL, 10);/* atoi replacement */
          free(uq_char);
          return uq;
     }          
}


void apply_uniq_threshold(var_t *var, uniq_filter_t *uniq_filter)
{
     if (! uniq_filter->thresh) {
          return;
     }

     if (uniq_phred_from_var(var) < uniq_filter->thresh) {
          vcf_var_add_to_filter(var, uniq_filter->id);
     }
}


/* returns -1 on error 
 *
 * filter everything that's not significant
 * 
 * FIXME should be part of lofreq filter.
 *
 */
int 
apply_uniq_filter_mtc(uniq_filter_t *uniq_filter, var_t **vars, const int num_vars)
{
     double *uniq_probs = NULL;
     int i;

     if (uniq_filter->ntests && num_vars > uniq_filter->ntests) {
         LOG_WARN("%s\n", "Number of predefined tests for uniq filter larger than number of variants! Are you sure that makes sense?");
     }

     if (! uniq_filter->ntests) {
          uniq_filter->ntests = num_vars;
     }

     /* collect uniq error probs
      */
     uniq_probs = malloc(num_vars * sizeof(double));
     if ( ! uniq_probs) {
          LOG_FATAL("%s\n", "out of memory");
          exit(1);
     }
     for (i=0; i<num_vars; i++) {
          uniq_probs[i] = PHREDQUAL_TO_PROB(uniq_phred_from_var(vars[i]));
     }

     /* multiple testing correction
      */
     if (uniq_filter->mtc_type == MTC_BONF) {
          bonf_corr(uniq_probs, num_vars, 
                    uniq_filter->ntests);
          
     } else if (uniq_filter->mtc_type == MTC_HOLMBONF) {
          holm_bonf_corr(uniq_probs, num_vars, 
                         uniq_filter->alpha, uniq_filter->ntests);
          
     } else if (uniq_filter->mtc_type == MTC_FDR) {
          int num_rej = 0;
          long int *idx_rej; /* indices of rejected i.e. significant values */
          int i;
          
          num_rej = fdr(uniq_probs, num_vars, 
                        uniq_filter->alpha, uniq_filter->ntests, 
                        &idx_rej);
          for (i=0; i<num_rej; i++) {
               int idx = idx_rej[i];
               uniq_probs[idx] = -1;
          }
          free(idx_rej);
          
     } else {
          LOG_FATAL("Internal error: unknown MTC type %d\n", uniq_filter->mtc_type);
          return -1;
     }

     for (i=0; i<num_vars; i++) {
          if (uniq_probs[i] > uniq_filter->alpha) {
               vcf_var_add_to_filter(vars[i], uniq_filter->id);
          }
     }

     free(uniq_probs);

     return 0;
}



/* used as pileup callback function which is not ideal since this can
 * only work on one position (has to be ensured by caller).
 *
 * No cov means I won't be called through mpileup and no output will
 * be generated. Non-sig pv means I'm not sure and no ouput will be
 * generated. Only if pv is sig we will print the var
 *
 * needs to return void to be used as function pointer to mpileup
 */
void
uniq_snv(const plp_col_t *p, void *confp)
{
     uniq_conf_t *conf = (uniq_conf_t *)confp;
     char *af_char = NULL;
     float af;
     int is_uniq = 0;
     int is_indel;
     int coverage;

     is_indel =  vcf_var_is_indel(conf->var);

#ifdef DISABLE_INDELS
     if (is_indel) {
          LOG_WARN("uniq logic can't be applied to indels."
                   " Skipping indel var at %s %d\n",
                   conf->var->chrom, conf->var->pos+1);
          return;
     }
#endif

     if (0 != strcmp(p->target, conf->var->chrom) || p->pos != conf->var->pos) {
          LOG_ERROR("wrong pileup for var. pileup for %s %d. var for %s %d\n",
                    p->target, p->pos+1, conf->var->chrom, conf->var->pos+1);
          return;
     }

     coverage = p->coverage_plp;
     if (is_indel) {
          coverage -= p->num_tails;
     }
     if (1 > coverage) {
          return;
     }

     if (conf->uni_freq <= 0.0) {
          if (! vcf_var_has_info_key(&af_char, conf->var, "AF")) {
               LOG_FATAL("%s\n", "Couldn't parse AF (key not found) from variant");
               /* hard to catch error later */
               exit(1);
          }
          af = strtof(af_char, (char **)NULL); /* atof */
          free(af_char);
          if (af < 0.0 || af > 1.0) {
               float new_af;
               new_af = af<0.0 ? 0.01 : 1.0;
               /* hard to catch error later */
               LOG_FATAL("Invalid (value out of bound) AF %f in variant. Resetting to %f\n", af, new_af);
               af = new_af;
          }

     } else {
          assert(conf->uni_freq <= 1.0);
          af = conf->uni_freq;
     }


     if (conf->use_det_lim) {
          /* given the current base counts and their error probs,
           * would we've been able to detect at given frequency.
           */
          long double pvalues[NUM_NONCONS_BASES];
          double *err_probs; /* error probs (qualities) passed down to snpcaller */
          int num_err_probs;

          int alt_bases[NUM_NONCONS_BASES];/* actual alt bases */
          int alt_counts[NUM_NONCONS_BASES]; /* counts for alt bases handed down to snpcaller */
          int alt_raw_counts[NUM_NONCONS_BASES]; /* raw, unfiltered alt-counts */
          varcall_conf_t varcall_conf;

          int bonf = 1;
          float alpha = 0.01;

          init_varcall_conf(&varcall_conf);
          if (debug) {
               dump_varcall_conf(&varcall_conf, stderr);
          }

          plp_to_errprobs(&err_probs, &num_err_probs,
                          alt_bases, alt_counts, alt_raw_counts,
                          p, &varcall_conf);
          LOG_DEBUG("at %s:%d with cov %d and num_err_probs %d\n", 
              p->target, p->pos, coverage, num_err_probs);

          /* Now pretend we see AF(SNV-to-test)*coverage variant
           * bases. Truncate to int, i.e err on the side of caution
           * during rounding (assume fewer alt bases) */
          alt_counts[0] = af * num_err_probs; /* don't use coverage as that is before filtering */
          alt_counts[1] = alt_counts[2] = 0;

          if (snpcaller(pvalues, err_probs, num_err_probs,
                        alt_counts, bonf, alpha)) {
               fprintf(stderr, "FATAL: snpcaller() failed at %s:%s():%d\n",
                       __FILE__, __FUNCTION__, __LINE__);
               free(err_probs);
               return;
          }

          /* only need to test first pv */
          if (pvalues[0] * (float)bonf < alpha) {
              /* significant value means given the counts and
               * qualities we would have been able to detect this
               * uncalled SNV had it been present at the given
               * frequency. But since we didn't this is a uniq
               * variant.
               * 
               * No point in adding this as phred qual because it
               * means the opposite of UQ
               */

               vcf_var_add_to_info(conf->var, uniq_flag);
          }

          LOG_VERBOSE("%s %d num_quals=%d assumed-var-counts=%d would-have-been-detectable=%d\n",
               conf->var->chrom, conf->var->pos+1, num_err_probs, alt_counts[0], is_uniq);
          free(err_probs);
          
     } else {
          int alt_count;
          double pvalue;
          char info_str[128];

          if (is_indel) {
               int ref_len = strlen(conf->var->ref);
               int alt_len = strlen(conf->var->alt);
               if (ref_len > alt_len) { /* deletion */
                    char *del_key = malloc((strlen(conf->var->ref)+1)*sizeof(char));
                    strcpy(del_key, conf->var->ref+1);
                    del_event *it_del = find_del_sequence(&p->del_event_counts, del_key);
                    if (it_del) {
                         alt_count = it_del->count;
                    } else {
                         alt_count = 0;
                    }
                    /* LOG_DEBUG("%s>%s k:%s c:%d\n", conf->var->ref, conf->var->alt, del_key, alt_count); */
                    free(del_key);
               } else { /* insertion */
                    char *ins_key = malloc((strlen(conf->var->alt)+1)*sizeof(char));
                    strcpy(ins_key, conf->var->alt+1);
                    ins_event *it_ins = find_ins_sequence(&p->ins_event_counts, ins_key);
                    if (it_ins) {
                         alt_count = it_ins->count;
                    } else {
                         alt_count = 0;
                    }
                    /* LOG_DEBUG("%s>%s k:%s c:%d\n", conf->var->ref, conf->var->alt, ins_key, alt_count);*/
                    free(ins_key);
               }

          } else {
               alt_count = base_count(p, conf->var->alt[0]);
          }


#ifdef DEBUG
          LOG_DEBUG("Now testing af=%f cov=%d alt_count=%d at %s %d for var:",
                    af, coverage, alt_count, p->target, p->pos+1);
#endif
          
          /* this is a one sided test */
          if (0 != binom(&pvalue, NULL, coverage, alt_count, af)) {
               LOG_ERROR("%s\n", "binom() failed");
               return;
          }

          snprintf(info_str, 128, "%s=%d", uniq_phred_tag, PROB_TO_PHREDQUAL_SAFE(pvalue));
          vcf_var_add_to_info(conf->var, info_str);

          LOG_DEBUG("%s %d %s>%s AF=%f | %s (p-value=%g) | BAM alt_count=%d cov=%d (freq=%f)\n",
                      conf->var->chrom, conf->var->pos+1, conf->var->ref, conf->var->alt, af,
                      is_uniq ? "unique" : "not necessarily unique", pvalue,
                      alt_count, coverage, alt_count/(float)coverage);
     }
}


static void
usage(const uniq_conf_t* uniq_conf)
{
     fprintf(stderr,
                  "\n%s: Checks whether variants predicted in one sample (listed in vcf input)" \
                  " are unique to this sample or if they were not called in other sample due" \
                  " to coverage issues. This is done by using a Binomial test with alternate"\
                  " and reference counts from the BAM and the variant frequency (i.e it's testing"\
                  " differences in frequencies. Alternatively, the logic can be changed to" \
                  " check whether the variant frequency would have been above LoFreq's" \
                  " detection limit given the BAM coverage and base-qualities."\
                  "\n\n" \
                  "Assigns UNIQ tag to variants considered unique."\
                  " Will ignore filtered input variants and will by default only report uniq variants.\n\n", MYNAME);

     fprintf(stderr,"Usage: %s [options] indexed-in.bam\n\n", MYNAME);
     fprintf(stderr,"Options:\n");
     fprintf(stderr, "  -v | --vcf-in FILE      Input vcf file listing variants [- = stdin; gzip supported]\n");
     fprintf(stderr, "  -o | --vcf-out FILE     Output vcf file [- = stdout; gzip supported]\n");
     fprintf(stderr, "  -f | --uni-freq         Assume variants have uniform test frequency of this value (unused if <=0) [%f]\n", uniq_conf->uni_freq);
     fprintf(stderr, "  -t | --uniq-thresh INT  Minimum uniq phred-value required. Conflicts with -m. 0 for off (default=%d)\n", uniq_conf->uniq_filter.thresh);
     fprintf(stderr, "  -m | --uniq-mtc STRING  Uniq multiple testing correction type. One of 'bonf', 'holm' or 'fdr'. (default=%s)\n", mtc_type_str[uniq_conf->uniq_filter.mtc_type]);
     fprintf(stderr, "  -a | --uniq-alpha FLOAT Uniq Multiple testing correction p-value threshold (default=%f)\n", uniq_conf->uniq_filter.alpha); 
     fprintf(stderr, "  -n | --uniq-ntests INT  Uniq multiple testing correction p-value threshold (default=#vars)\n");
     fprintf(stderr, "       --output-all       Report all variants instead of only the ones, marked unique.\n");
     fprintf(stderr, "                          Note, that variants already filtered in input will not be printed.\n");
     fprintf(stderr, "       --use-det-lim      Report variants if they are above implied detection limit\n");
     fprintf(stderr, "                          Default is to use binomial test to check for frequency differences\n");
     fprintf(stderr, "       --use-orphan       Don't ignore anomalous read pairs / orphan reads\n");
     fprintf(stderr, "       --verbose          Be verbose\n");
     fprintf(stderr, "       --debug            Enable debugging\n");
}
/* usage() */


int
main_uniq(int argc, char *argv[])
{
     int c, i;
     char *bam_file = NULL;
     char *vcf_in = NULL; /* - == stdout */
     char *vcf_out = NULL; /* - == stdout */
     mplp_conf_t mplp_conf;
     uniq_conf_t uniq_conf;
     void (*plp_proc_func)(const plp_col_t*, void*);
     int rc = 0;
     var_t **vars = NULL;
     int num_vars = 0;
     char *vcf_header = NULL;
     static int use_det_lim = 0;
     static int use_orphan = 0;
     static int output_all = 0;
     static int is_somatic = 0;

     /* default uniq options */
     memset(&uniq_conf, 0, sizeof(uniq_conf_t));
     uniq_conf.uni_freq = DEFAULT_UNI_FREQ;
     uniq_conf.use_det_lim = 0;

     uniq_conf.uniq_filter.mtc_type = MTC_FDR;
     uniq_conf.uniq_filter.alpha = 0.001;

     /* default pileup options */
     memset(&mplp_conf, 0, sizeof(mplp_conf_t));
     mplp_conf.max_mq = DEFAULT_MAX_MQ;
     mplp_conf.min_mq = 1;
     mplp_conf.min_plp_bq = DEFAULT_MIN_PLP_BQ;
     mplp_conf.max_depth = DEFAULT_MAX_PLP_DEPTH;
     mplp_conf.flag = MPLP_NO_ORPHAN;


    /* keep in sync with long_opts_str and usage
     *
     * getopt is a pain in the whole when it comes to syncing of long
     * and short args and usage. check out gopt, libcfu...
     */
    while (1) {
         static struct option long_opts[] = {
              /* see usage sync */
              {"help", no_argument, NULL, 'h'},
              {"verbose", no_argument, &verbose, 1},
              {"debug", no_argument, &debug, 1},
              {"use-det-lim", no_argument, &use_det_lim, 1},
              {"use-orphan", no_argument, &use_orphan, 1},
              {"output-all", no_argument, &output_all, 1},
              {"is-somatic", no_argument, &is_somatic, 1},

              {"vcf-in", required_argument, NULL, 'v'},
              {"vcf-out", required_argument, NULL, 'o'},

              {"uni-freq", required_argument, NULL, 'f'},

              {"uniq-thresh", required_argument, NULL, 't'},
              {"uniq-mtc", required_argument, NULL, 'm'},
              {"uniq-alpha", required_argument, NULL, 'a'},
              {"uniq-ntests", required_argument, NULL, 'n'},

              {0, 0, 0, 0} /* sentinel */
         };

         /* keep in sync with long_opts and usage */
         static const char *long_opts_str = "hv:o:f:t:m:a:n:";

         /* getopt_long stores the option index here. */
         int long_opts_index = 0;
         c = getopt_long(argc-1, argv+1, /* skipping 'lofreq', just leaving 'command', i.e. call */
                         long_opts_str, long_opts, & long_opts_index);
         if (c == -1) {
              break;
         }

         switch (c) {
         /* keep in sync with long_opts etc */
         case 'h':
              usage(& uniq_conf);
              return 0;

         case 'v':
              if (0 != strcmp(optarg, "-")) {
                   if (! file_exists(optarg)) {
                        LOG_FATAL("Input file '%s' does not exist. Exiting...\n", optarg);
                        return 1;
                   }
              }
              vcf_in = strdup(optarg);
              break;

         case 'o':
              if (0 != strcmp(optarg, "-")) {
                   if (file_exists(optarg)) {
                        LOG_FATAL("Cowardly refusing to overwrite file '%s'. Exiting...\n", optarg);
                        return 1;
                   }
              }
              vcf_out = strdup(optarg);
              break;

         case 'f':
              uniq_conf.uni_freq = strtof(optarg, (char **)NULL); /* atof */
              if (uniq_conf.uni_freq<=0) {
                   LOG_WARN("%s\n", "Ignoring uni-freq option");
              }
              if (uniq_conf.uni_freq>1.0) {
                   LOG_FATAL("%s\n", "Value for uni-freq has to be <1.0");
                   return 1;
              }
              break;

         case 't':
              if (! isdigit(optarg[0])) {
                   LOG_FATAL("Non-numeric argument provided: %s\n", optarg);
                   return -1;
              }
              uniq_conf.uniq_filter.thresh = atoi(optarg);
              uniq_conf.uniq_filter.mtc_type = MTC_NONE;
              break;
         case 'm':
              uniq_conf.uniq_filter.mtc_type = mtc_str_to_type(optarg);
              if (-1 == uniq_conf.uniq_filter.mtc_type) {
                   LOG_FATAL("Unknown multiple testing correction type '%s' for snv quality filtering\n", optarg);
                   return -1;
              }
              break;
         case 'a':
              if (! isdigit(optarg[0])) {
                   LOG_FATAL("Non-numeric argument provided: %s\n", optarg);
                   return -1;
              }
              uniq_conf.uniq_filter.alpha = strtof(optarg, NULL);
              break;
         case 'n':
              if (! isdigit(optarg[0])) {
                   LOG_FATAL("Non-numeric argument provided: %s\n", optarg);
                   return -1;
              }
              uniq_conf.uniq_filter.ntests = atol(optarg);
              break;

         case '?':
              LOG_FATAL("%s\n", "unrecognized arguments found. Exiting...\n");
              return 1;
         default:
              break;
         }
    }
    if (use_orphan) {
         mplp_conf.flag &= ~MPLP_NO_ORPHAN;
    }
    if (debug) {
         dump_mplp_conf(& mplp_conf, stderr);
    }
    uniq_conf.output_all = output_all;
    uniq_conf.use_det_lim = use_det_lim;


#if DEBUG
    LOG_DEBUG("uniq_conf.uniq_filter.thresh = %d\n", uniq_conf.uniq_filter.thresh);
    LOG_DEBUG("uniq_conf.uniq_filter.mtc_type = %d\n", uniq_conf.uniq_filter.mtc_type);
    LOG_DEBUG("uniq_conf.uniq_filter.alpha = %f\n", uniq_conf.uniq_filter.alpha);
    LOG_DEBUG("uniq_conf.uniq_filter.ntests = %d\n", uniq_conf.uniq_filter.ntests);
#endif
    
    if (uniq_conf.uniq_filter.thresh && uniq_conf.uniq_filter.mtc_type != MTC_NONE) {
         LOG_FATAL("%s\n", "Can't use fixed Unique quality threshold *and* multiple testing correction.");
         return 1;
    }

    if (argc == 2) {
        fprintf(stderr, "\n");
        usage(& uniq_conf);
        return 1;
    }

    if (1 != argc - optind - 1) {
        fprintf(stderr, "Need exactly one BAM file as last argument\n");
        return 1;
    }
    bam_file = (argv + optind + 1)[0];
    if (! file_exists(bam_file)) {
         LOG_FATAL("BAM file %s does not exist. Exiting...\n", bam_file);
         return -1;
    }


    if (! vcf_in) {
#if 0
         vcf_in = malloc(2 * sizeof(char));
         strcpy(vcf_in, "-");
#else
         LOG_FATAL("%s\n", "No input vcf specified. Exiting...");
         return -1;
#endif
    }
    if (! vcf_out) {
         vcf_out = malloc(2 * sizeof(char));
         strcpy(vcf_out, "-");
    }

    if (vcf_file_open(& uniq_conf.vcf_in, vcf_in,
                      HAS_GZIP_EXT(vcf_in), 'r')) {
         LOG_ERROR("Couldn't open %s\n", vcf_in);
         return 1;
    }

    if (vcf_file_open(& uniq_conf.vcf_out, vcf_out,
                      HAS_GZIP_EXT(vcf_out), 'w')) {
         LOG_ERROR("Couldn't open %s\n", vcf_out);
         return 1;
    }

    if (0 != vcf_parse_header(&vcf_header, & uniq_conf.vcf_in)) {
         LOG_WARN("%s\n", "vcf_parse_header() failed. trying to rewind to start...");
         if (vcf_file_seek(& uniq_conf.vcf_in, 0, SEEK_SET)) {
              LOG_FATAL("%s\n", "Couldn't rewind file to parse variants"
                        " after header parsing failed");
              return 1;
         }
    } else {
         vcf_header_add(&vcf_header, "##INFO=<ID=UNIQ,Number=0,Type=Flag,Description=\"Unique, i.e. not detectable in paired sample\">\n");
         vcf_header_add(&vcf_header, "##INFO=<ID=UQ,Number=1,Type=Integer,Description=\"Phred-scaled uniq score at this position\">\n");
         if (is_somatic) {
              vcf_header_add(&vcf_header, "##INFO=<ID=SOMATIC,Number=0,Type=Flag,Description=\"Somatic event\">\n");
         }
         if (! uniq_conf.use_det_lim) {
              char full_filter_str[FILTER_STRSIZE];
              if (uniq_conf.uniq_filter.thresh > 0) {
                   snprintf(uniq_conf.uniq_filter.id, FILTER_ID_STRSIZE, "min_uq_%d", uniq_conf.uniq_filter.thresh);
                   snprintf(full_filter_str, FILTER_STRSIZE,
                            "##FILTER=<ID=%s,Description=\"Minimum Uniq Phred %d\">\n",
                            uniq_conf.uniq_filter.id, uniq_conf.uniq_filter.thresh);
                   vcf_header_add(&vcf_header, full_filter_str);
                   
              } else if (uniq_conf.uniq_filter.mtc_type != MTC_NONE) {
                   char buf[64];
                   mtc_str(buf, uniq_conf.uniq_filter.mtc_type);
                   snprintf(uniq_conf.uniq_filter.id, FILTER_ID_STRSIZE, "uq_%s", buf);
                   snprintf(full_filter_str, FILTER_STRSIZE,
                            "##FILTER=<ID=%s,Description=\"Uniq Multiple Testing Correction: %s corr. pvalue < %f\">\n",
                            uniq_conf.uniq_filter.id, buf, uniq_conf.uniq_filter.alpha);
                   vcf_header_add(& vcf_header, full_filter_str);
              }
         }

         vcf_write_header(& uniq_conf.vcf_out, vcf_header);
         free(vcf_header);
    }

    num_vars = vcf_parse_vars(&vars, & uniq_conf.vcf_in, 1);
    if (0 == num_vars) {
         LOG_WARN("%s\n", "Didn't find any variants in input");
         goto clean_and_exit;
    }
    if (! uniq_conf.uniq_filter.ntests) {
         uniq_conf.uniq_filter.ntests = num_vars;
    }

    plp_proc_func = &uniq_snv;

    for (i=0; i<num_vars; i++) {
         char reg_buf[BUF_SIZE];
         if (i%100==0) {
              LOG_VERBOSE("Processing variant %d of %d\n", i+1, num_vars);
         }
         uniq_conf.var = vars[i];

         snprintf(reg_buf, BUF_SIZE, "%s:%ld-%ld",
                  vars[i]->chrom, vars[i]->pos+1, vars[i]->pos+1);
         mplp_conf.reg = strdup(reg_buf);

         LOG_DEBUG("pileup for var no %d at %s %d\n",
                   i+1, uniq_conf.var->chrom, uniq_conf.var->pos+1);
#ifdef DISABLE_INDELS
         if (vcf_var_has_info_key(NULL, uniq_conf.var, "INDEL")) {
              LOG_WARN("Skipping indel var at %s %d\n",
                       uniq_conf.var->chrom, uniq_conf.var->pos+1);
              free(mplp_conf.reg);
              mplp_conf.reg = NULL;
              continue;
         }
#endif
         /* no need to check for filter because done by parse_vars */

         rc = mpileup(&mplp_conf, plp_proc_func, (void*)&uniq_conf,
                      1, (const char **) argv + optind + 1);

         if (uniq_conf.uniq_filter.thresh) {
              apply_uniq_threshold(uniq_conf.var, & uniq_conf.uniq_filter);
         }

         free(mplp_conf.reg);
         mplp_conf.reg = NULL;
    }
    uniq_conf.var = NULL;/* just be sure to not use it accidentally again */


    /* print whatever we've got. there's no UQ to test or we
     * are supposed to print all 
     */
    if (uniq_conf.use_det_lim) {
         for (i=0; i<num_vars; i++) {
              var_t *var = vars[i];
              vcf_write_var(& uniq_conf.vcf_out, var);
         }
         /* all done */
         goto clean_and_exit;
    }



    if (uniq_conf.uniq_filter.mtc_type != MTC_NONE) {
         if (apply_uniq_filter_mtc(& uniq_conf.uniq_filter, vars, num_vars)) {
              LOG_FATAL("%s\n", "Multiple testing correction on uniq pvalues failed");
              return -1;
         }
    }
    
    for (i=0; i<num_vars; i++) {
         var_t *var = vars[i];
         if (VCF_VAR_PASSES(var) || uniq_conf.output_all) {
              vcf_write_var(& uniq_conf.vcf_out, var);
         }
    }

clean_and_exit:

    vcf_file_close(& uniq_conf.vcf_in);
    vcf_file_close(& uniq_conf.vcf_out);

    for (i=0; i<num_vars; i++) {
         vcf_free_var(& vars[i]);
    }
    free(vars);

    free(vcf_in);
    free(vcf_out);

    if (0==rc) {
         LOG_VERBOSE("%s\n", "Successful exit.");
    }
    /* LOG_FIXME("%s\n", "allow user setting of -S and -J. Currently just using default") */

    return rc;
}
/* main_uniq */

