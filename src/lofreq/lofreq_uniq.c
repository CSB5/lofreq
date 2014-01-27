/* -*- c-file-style: "k&r"; indent-tabs-mode: nil; -*- */


/*********************************************************************
 *
 * FIXME update license
 *
 *********************************************************************/


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


#if 1
#define MYNAME "lofreq uniq"
#else
#define MYNAME PACKAGE
#endif

#define DEFAULT_SIG 0.05

#define DEFAULT_UNI_FREQ -1.0

#define BUF_SIZE 1<<16

typedef struct {
     /* fix values once established */
     float sig;
     float uni_freq;
     vcf_file_t vcf_out;
     vcf_file_t vcf_in;
     long long int bonf;
     int use_freq_diff;

     /* changing per pos: the var to test */
     var_t *var;
} uniq_conf_t;



/* used as pileup callback function which is not ideal since this can
 * only work on one position (has to be ensured by caller).
 *
 * No cov means I won't be called and no output will be generated.
 * Non-sig pv means I'm not sure and no ouput will be generated. Only
 * if pv is sig we will print the var
 */
void
uniq_snv(const plp_col_t *p, void *confp)
{
     uniq_conf_t *conf = (uniq_conf_t *)confp;
     char *af_char = NULL;
     float af;
     int is_uniq = 0;

     if (vcf_var_has_info_key(NULL, conf->var, "INDEL")) {
          LOG_WARN("uniq logic can't be applied to indels."
                   " Skipping indel var at %s %d\n",
                   conf->var->chrom, conf->var->pos+1);
          return;
     }

     if (0 != strcmp(p->target, conf->var->chrom) || p->pos != conf->var->pos) {
          LOG_ERROR("wrong pileup for var. pileup for %s %d. var for %s %d\n",
                    p->target, p->pos+1, conf->var->chrom, conf->var->pos+1);
          return;
     }

     /* no coverage usually means I won't be called, but to be safe: */
     if (0 == p->coverage) {
          return;
     }

     if (conf->uni_freq <= 0.0) {
          vcf_var_has_info_key(&af_char, conf->var, "AF");
          if (NULL == af_char) {
#if 0
               LOG_ERROR("%s\n", "Couldn't parse AF (key not found) from the following variant:");
               vcf_write_var(stderr, conf->var);
#else
               LOG_ERROR("%s\n", "Couldn't parse AF (key not found) from variant");
#endif
               return;
          }
          af = strtof(af_char, (char **)NULL); /* atof */
          free(af_char);
          if (af < 0.0 || af > 1.0) {
#if 0
               LOG_ERROR("%s\n", "Couldn't parse AF (value out of bound) from the following variant:");
               vcf_write_var(stderr, conf->var);
#else
               LOG_ERROR("%s\n", "Couldn't parse AF (value out of bound) from variant");
#endif
               return;
          }
     } else {
          assert(conf->uni_freq <= 1.0);
          af = conf->uni_freq;
     }


    if (! conf->use_freq_diff) {
          /* given the current base counts and their error probs,
           * would we've been able to detect at given frequency.
           */
          double pvalues[NUM_NONCONS_BASES];
          double *err_probs; /* error probs (qualities) passed down to snpcaller */
          int num_err_probs;

          int alt_bases[NUM_NONCONS_BASES];/* actual alt bases */
          int alt_counts[NUM_NONCONS_BASES]; /* counts for alt bases handed down to snpcaller */
          int alt_raw_counts[NUM_NONCONS_BASES]; /* raw, unfiltered alt-counts */
          snvcall_conf_t snvcall_conf;

          init_snvcall_conf(&snvcall_conf);
          if (debug) {
               dump_snvcall_conf(&snvcall_conf, stderr);
          }

          plp_to_errprobs(&err_probs, &num_err_probs,
                          alt_bases, alt_counts, alt_raw_counts,
                          p, &snvcall_conf);
          LOG_DEBUG("at %s:%d with cov %d and num_err_probs %d\n", 
              p->target, p->pos, p->coverage, num_err_probs);

          /* Now pretend we see AF(SNV-to-test)*coverage variant
           * bases. Truncate to int, i.e err on the side of caution
           * during rounding (assume fewer alt bases) */
          alt_counts[0] = af * num_err_probs; /* don't use p->coverage as that is before filtering */
          alt_counts[1] = alt_counts[2] = 0;

          if (snpcaller(pvalues, err_probs, num_err_probs,
                        alt_counts, conf->bonf, conf->sig)) {
               fprintf(stderr, "FATAL: snpcaller() failed at %s:%s():%d\n",
                       __FILE__, __FUNCTION__, __LINE__);
               free(err_probs);
               return;
          }

          /* only need to test first one */
          if (pvalues[0] * (float)conf->bonf < conf->sig) {
              /* sign. value, i.e. we would have been able to detect
               * this uncalled SNV, but didn't */
              is_uniq = 1;
          }

          LOG_VERBOSE("%s %d num_quals=%d assumed-var-counts=%d would-have-been-detectable=%d\n",
               conf->var->chrom, conf->var->pos+1, num_err_probs, alt_counts[0], is_uniq);
          free(err_probs);
          
     } else {
          int alt_count = base_count(p, conf->var->alt);
          double pvalue;

#ifdef DEBUG
          LOG_DEBUG("Now testing af=%f cov=%d alt_count=%d at %s %d for var:",
          af, p->coverage, alt_count, p->target, p->pos+1);
#endif
          
          /* this is a one sided test */
          if (0 != binom(&pvalue, NULL, p->coverage, alt_count, af)) {
               LOG_ERROR("%s\n", "binom() failed");
               return;
          }

          if (pvalue < conf->sig/(float)conf->bonf) {
               char info_str[128];
               is_uniq = 1;
               snprintf(info_str, 128, "UQ=%d", PROB_TO_PHREDQUAL(pvalue));
               vcf_var_add_to_info(conf->var, info_str);
          }
          LOG_VERBOSE("%s %d %c>%c AF=%f | %s (p-value=%g sig/bonf=%g) | BAM alt_count=%d cov=%d (freq=%f)\n",
                      conf->var->chrom, conf->var->pos+1, conf->var->ref, conf->var->alt, af,
                      is_uniq ? "unique" : "not necessarily unique", pvalue, conf->sig/(float)conf->bonf,
                      alt_count, p->coverage, alt_count/(float)p->coverage);

     }

    if (is_uniq) {
        vcf_write_var(& conf->vcf_out, conf->var);
    }
}


static void
usage(const uniq_conf_t* uniq_conf)
{
     fprintf(stderr,
                  "%s: Checks whether variants predicted in one sample (listed in vcf input)" \
                  " are unique to this sample or if they were not called in other sample due" \
                  " to coverage issues. This is either done by" \
                  " (1) checking whether the variant frequency would have been above LoFreq's" \
                  " detection limit given the BAM coverage and base-qualities or" \
                  " (2) using a Binomial the BAM alternate and reference counts and the variant frequency (more suited for testing freq differences)." \
                  "\n\n" \
                  "Will ignore filtered input variants as well as indels and will only" \
                  " output variants considered unique.\n\n", MYNAME);

     fprintf(stderr,"Usage: %s [options] indexed-in.bam\n\n", MYNAME);
     fprintf(stderr,"Options:\n");
     fprintf(stderr, "  -v | --vcf-in FILE       Input vcf file listing variants [- = stdin; gzip supported]\n");
     fprintf(stderr, "  -o | --vcf-out FILE      Output vcf file [- = stdout; gzip supported]\n");
     fprintf(stderr, "  -s | --sig               Significance threshold [%f]\n", uniq_conf->sig);
     fprintf(stderr, "  -f | --uni-freq          Assume variants have uniform test frequency of this value (unused if <=0) [%f]\n", uniq_conf->uni_freq);
     fprintf(stderr, "       --use-freq-diff     Use binomial test to check for frequency differences\n");
     fprintf(stderr, "                           (Probably not a good idea at high coverage positions)\n");
     fprintf(stderr, "                           Default is to report variants if they are above implied detection limit\n");
     fprintf(stderr, "       --use-orphan        Count anomalous read pairs\n");
     fprintf(stderr, "       --verbose           Be verbose\n");
     fprintf(stderr, "       --debug             Enable debugging\n");
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
     static int use_freq_diff = 0;
     static int use_orphan = 0;

     for (i=0; i<argc; i++) {
          LOG_DEBUG("arg %d: %s\n", i, argv[i]);
     }


     /* default uniq options */
     memset(&uniq_conf, 0, sizeof(uniq_conf_t));
     uniq_conf.sig = DEFAULT_SIG;
     uniq_conf.uni_freq = DEFAULT_UNI_FREQ;
     uniq_conf.bonf = 1;
     uniq_conf.use_freq_diff = 0;

     /* default pileup options */
     memset(&mplp_conf, 0, sizeof(mplp_conf_t));
     mplp_conf.max_mq = DEFAULT_MAX_MQ;
     mplp_conf.min_mq = 1;
     mplp_conf.min_bq = DEFAULT_MIN_BQ;
     mplp_conf.capQ_thres = 0;
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
              {"use-freq-diff", no_argument, &use_freq_diff, 1},
              {"use-orphan", no_argument, &use_orphan, 1},

              {"vcf-in", required_argument, NULL, 'v'},
              {"vcf-out", required_argument, NULL, 'o'},
              {"uni-freq", required_argument, NULL, 'f'},
              {"sig", required_argument, NULL, 's'},

              {0, 0, 0, 0} /* sentinel */
         };

         /* keep in sync with long_opts and usage */
         static const char *long_opts_str = "hv:o:s:f:";

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

         case 's':
              uniq_conf.sig = strtof(optarg, (char **)NULL); /* atof */
              if (0==uniq_conf.sig) {
                   LOG_FATAL("%s\n", "Couldn't parse sign-threshold");
                   return 1;
              }
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

    uniq_conf.use_freq_diff = use_freq_diff;

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
              return -1;
         }
    } else {
         /*LOG_FIXME("header was: %s", vcf_header);*/
         vcf_header_add(&vcf_header, "##INFO=<ID=UQ,Number=1,Type=Integer,Description=\"Phred-scaled uniq score at this position\">\n");
         vcf_write_header(& uniq_conf.vcf_out, vcf_header);
         /*LOG_FIXME("header now: %s", vcf_header);*/
         free(vcf_header);
    }


    if (-1 == (num_vars = vcf_parse_vars(&vars, & uniq_conf.vcf_in, 1))) {
         LOG_FATAL("%s\n", "vcf_parse_vars() failed");
         return -1;
    }

    plp_proc_func = &uniq_snv;
    uniq_conf.bonf = num_vars;
    if (0 == num_vars) {
         LOG_WARN("%s\n", "Didn't find any variants in input");
    }
    for (i=0; i<num_vars; i++) {
         char reg_buf[BUF_SIZE];

         uniq_conf.var = vars[i];

         snprintf(reg_buf, BUF_SIZE, "%s:%ld-%ld",
                  vars[i]->chrom, vars[i]->pos+1, vars[i]->pos+1);
         mplp_conf.reg = strdup(reg_buf);

         LOG_DEBUG("pileup for var no %d at %s %d\n",
                   i+1, uniq_conf.var->chrom, uniq_conf.var->pos+1);
         if (vcf_var_has_info_key(NULL, uniq_conf.var, "INDEL")) {
              LOG_VERBOSE("Skipping indel var at %s %d\n",
                       uniq_conf.var->chrom, uniq_conf.var->pos+1);
              free(mplp_conf.reg);
              mplp_conf.reg = NULL;
              continue;

         } else if (vcf_var_filtered(uniq_conf.var)) {
              LOG_VERBOSE("Skipping filtered var at %s %d\n",
                       uniq_conf.var->chrom, uniq_conf.var->pos+1);
              free(mplp_conf.reg);
              mplp_conf.reg = NULL;
              continue;
         }

         rc = mpileup(&mplp_conf, plp_proc_func, (void*)&uniq_conf,
                      1, (const char **) argv + optind + 1);

         free(mplp_conf.reg);
         mplp_conf.reg = NULL;
    }

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
/* main_call */
