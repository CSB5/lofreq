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


#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <errno.h>
#include <assert.h>
#include <stdarg.h>
#include <getopt.h>
#include <stdlib.h>

/* lofreq includes */
#include "lofreq_filter.h"
#include "vcf.h"
#include "log.h"
#include "utils.h"
#include "multtest.h"
#include "defaults.h"


#if 1
#define MYNAME "lofreq filter"
#else
#define MYNAME PACKAGE
#endif

#define FILTER_ID_STRSIZE 64
#define FILTER_STRSIZE 128

#define ALT_STRAND_RATIO 0.85

typedef struct {
     int min;
     char id_min[FILTER_ID_STRSIZE];
     int max;
     char id_max[FILTER_ID_STRSIZE];
} dp_filter_t;

typedef struct {
     float min;
     char id_min[FILTER_ID_STRSIZE];
     float max;
     char id_max[FILTER_ID_STRSIZE];
} af_filter_t;

typedef struct {
     int thresh;/* use if > 0; otherwise use multiple testing correction that's if >0 */
     int mtc_type;/* holm; holmbonf; fdr; none */
     double alpha;
     long int ntests;
     char id[FILTER_ID_STRSIZE];
     int no_compound; /* otherwise ALT_STRAND_RATIO of var bases have to be on one strand as well */
     int incl_indels; /* if 1, also apply to indels */
} sb_filter_t;

typedef struct {
     int thresh;/* use if > 0; otherwise use multiple testing correction that's if >0 */
     int mtc_type;/* holm; holmbonf; fdr; none */
     double alpha;
     long int ntests;
     char id[FILTER_ID_STRSIZE];
} snvqual_filter_t;

typedef struct {
     int thresh;/* use if > 0; otherwise use multiple testing correction that's if >0 */
     int mtc_type;/* holm; holmbonf; fdr; none */
     double alpha;
     long int ntests;
     char id[FILTER_ID_STRSIZE];
} indelqual_filter_t;

typedef struct {
     vcf_file_t vcf_in;
     vcf_file_t vcf_out;
     int print_only_passed;
     int only_snvs;
     int only_indels;

     /* each allowed to be NULL if not set */
     dp_filter_t dp_filter;
     af_filter_t af_filter;
     sb_filter_t sb_filter;
     snvqual_filter_t snvqual_filter;
     indelqual_filter_t indelqual_filter;
} filter_conf_t;


static int af_missing_warning_printed = 0;
static int dp_missing_warning_printed = 0;
static int dp4_missing_warning_printed = 0;
static int sb_missing_warning_printed = 0;


void
dump_filter_conf(const filter_conf_t *cfg)
{
     fprintf(stderr, "filter_conf:\n");
     fprintf(stderr, "  print_only_passed=%d\n", cfg->print_only_passed);
     fprintf(stderr, "  only_snvs=%d\n", cfg->only_snvs);
     fprintf(stderr, "  only_indels=%d\n", cfg->only_indels);

     fprintf(stderr, "  dp_filter min=%d max=%d\n",
             cfg->dp_filter.min, cfg->dp_filter.max);
     fprintf(stderr, "  af_filter min=%f max=%f\n",
             cfg->af_filter.min, cfg->af_filter.max);
     fprintf(stderr, "  sb_filter thresh=%d mtc_type=%d|%s alpha=%f ntests=%ld no_compound=%d incl_indel=%d\n",
             cfg->sb_filter.thresh, cfg->sb_filter.mtc_type, mtc_type_str[cfg->sb_filter.mtc_type],
             cfg->sb_filter.alpha, cfg->sb_filter.ntests, cfg->sb_filter.no_compound, cfg->sb_filter.incl_indels);
     fprintf(stderr, "  snvqual_filter thresh=%d mtc_type=%d|%s alpha=%f ntests=%ld\n",
             cfg->snvqual_filter.thresh, cfg->snvqual_filter.mtc_type, mtc_type_str[cfg->snvqual_filter.mtc_type],
             cfg->snvqual_filter.alpha, cfg->snvqual_filter.ntests);
     fprintf(stderr, "  indelqual_filter thresh=%d mtc_type=%d|%s alpha=%f ntests=%ld\n",
             cfg->indelqual_filter.thresh, cfg->indelqual_filter.mtc_type, mtc_type_str[cfg->indelqual_filter.mtc_type],
             cfg->indelqual_filter.alpha, cfg->indelqual_filter.ntests);
}


static void
usage(const filter_conf_t* filter_conf)
{
     fprintf(stderr, "%s: Filter variant parsed from vcf file\n\n", MYNAME);
     fprintf(stderr, "Usage: %s [options] -i input.vcf -o output.vcf\n", MYNAME);

     fprintf(stderr,"Options:\n");
     fprintf(stderr, "  Files:\n");
     fprintf(stderr, "  -i | --in FILE                 VCF input file (gzip supported)\n");
     fprintf(stderr, "  -o | --out FILE                VCF output file (default: - for stdout; gzip supported).\n");

     fprintf(stderr, "  Coverage (DP):\n");
     fprintf(stderr, "  -v | --cov-min INT             Minimum coverage allowed (<1=off)\n");
     fprintf(stderr, "  -V | --cov-max INT             Maximum coverage allowed (<1=off)\n");

     fprintf(stderr, "  Allele Frequency (AF; neg. values = off):\n");
     fprintf(stderr, "  -a | --af-min FLOAT            Minimum allele freq allowed (<1=off)\n");
     fprintf(stderr, "  -A | --af-max FLOAT            Maximum allele freq allowed (<1=off)\n");

     fprintf(stderr, "\n");
     fprintf(stderr, "  Strand Bias (SB):\n");
     fprintf(stderr, "  Note, variants are only filtered if their SB pvalue is below the threshold\n");
     fprintf(stderr, "  AND %d%% of variant bases are on one strand (toggled with --sb-no-compound).\n", (int)(ALT_STRAND_RATIO*100));
     fprintf(stderr, "  -B | --sb-thresh INT           Maximum phred-value allowed. Conflicts with -b.\n");
     fprintf(stderr, "  -b | --sb-mtc STRING           Multiple testing correction type. One of 'bonf', 'holm' or 'fdr'. Conflicts with -B\n");
     fprintf(stderr, "  -c | --sb-alpha FLOAT          Multiple testing correcion pvalue threshold\n");
     fprintf(stderr, "       --sb-no-compound          Don't use compound filter\n");
     fprintf(stderr, "       --sb-incl-indels          Apply SB filter to indels as well\n");

     fprintf(stderr, "\n");
     fprintf(stderr, "  SNV Quality:\n");
     fprintf(stderr, "  -Q | --snvqual-thresh INT      Minimum phred-value allowed. Conflicts with -q\n");
     fprintf(stderr, "  -q | --snvqual-mtc STRING      Multiple testing correction type. One of 'bonf', 'holm' or 'fdr'. Conflicts with -Q\n");
     fprintf(stderr, "  -r | --snvqual-alpha FLOAT     Multiple testing correcion pvalue threshold\n");
     fprintf(stderr, "  -s | --snvqual-ntests INT      Multiple testing correcion pvalue threshold\n");

     fprintf(stderr, "\n");
     fprintf(stderr, "  Indels:\n");
     fprintf(stderr, "  -K | --indelqual-thresh INT    Minimum phred-value allowed. Conflicts with -q\n");
     fprintf(stderr, "  -k | --indelqual-mtc STRING    Multiple testing correction type. One of 'bonf', 'holm' or 'fdr'. Conflicts with -Q\n");
     fprintf(stderr, "  -l | --indelqual-alpha FLOAT   Multiple testing correcion pvalue threshold\n");
     fprintf(stderr, "  -m | --indelqual-ntests INT    Multiple testing correcion pvalue threshold\n");

     fprintf(stderr, "\n");
     fprintf(stderr, "  Misc.:\n");
     fprintf(stderr, "       --only-indels             Keep InDels only\n");
     fprintf(stderr, "       --only-snvs               Keep SNVs only\n");
     fprintf(stderr, "       --print-all               Print all, not just passed variants\n");
     fprintf(stderr, "       --no-defaults             Remove all default filter settings\n");
     fprintf(stderr, "       --verbose                 Be verbose\n");
     fprintf(stderr, "       --debug                   Enable debugging\n");
     fprintf(stderr, "\nNOTE: without --no-defaults LoFreq's predefined filters are on (run with --verbose to see details)\n");
     fprintf(stderr, "\n");
}
/* usage() */



int alt_mostly_on_one_strand(var_t *var)
{
     dp4_counts_t dp4;
     float ratio = 0.0;

     if (vcf_get_dp4(&dp4, var)) {
          if (! dp4_missing_warning_printed) {
               LOG_WARN("%s\n", "DP4 info missing. Compound SB filter won't work");
               dp4_missing_warning_printed = 1;
          }
          return 0;
     }          

     /* FIXME: also check whether ref and alt ration is opposite?
        pro: that's the FPs we usually see
        con: violating fisher's exact test and additional rather arbitrary filter  */

     ratio = MAX(dp4.alt_fw, dp4.alt_rv)/(float)(dp4.alt_fw + dp4.alt_rv);
#if 0
     LOG_DEBUG("ratio for %s %d = %f\n", var->chrom, var->pos, ratio);
#endif
     if (ratio > ALT_STRAND_RATIO) {
          return 1;
     } else {
          return 0;
     }
}


void apply_af_filter(var_t *var, af_filter_t *af_filter)
{
     char *af_char = NULL;
     float af;

     if (af_missing_warning_printed) {
          return;
     }

     if (af_filter->min > 0 || af_filter->max > 0) {
          if ( ! vcf_var_has_info_key(&af_char, var, "AF")) {
               if ( ! af_missing_warning_printed) {
                    LOG_WARN("%s\n", "Requested AF filtering failed since AF tag is missing in variant");
                    af_missing_warning_printed = 1;
                    return;
               }
          }
          af = strtof(af_char, (char **)NULL); /* atof */
          if (errno==ERANGE) {
               LOG_ERROR("Couldn't parse EF from af_char %s. Disabling AF filtering", af_char);
               af_missing_warning_printed = 1;
               return;
          }
          free(af_char);

          if (af_filter->min > 0.0 && af < af_filter->min) {
               vcf_var_add_to_filter(var, af_filter->id_min);
          }
          if (af_filter->max > 0.0 && af > af_filter->max) {
               vcf_var_add_to_filter(var, af_filter->id_max);
          }
     }
}


void apply_dp_filter(var_t *var, dp_filter_t *dp_filter)
{
     char *dp_char = NULL;
     int cov;

     if (dp_missing_warning_printed) {
          return;
     }

     if (dp_filter->min > 0 || dp_filter->max > 0) {
          if ( ! vcf_var_has_info_key(&dp_char, var, "DP")) {
               if ( ! dp_missing_warning_printed) {
#ifdef DEBUG
                    vcf_file_t f; f.fh = stderr; f.gz = 0; vcf_write_var(&f, var);
#endif
                    LOG_WARN("%s\n", "Requested coverage filtering failed since DP tag is missing in variant");
                    dp_missing_warning_printed = 1;
                    return;
               }
          }
          errno = 0;
          /*cov = atoi(dp_char);*/
          cov = strtol(dp_char, (char **) NULL, 10);
          if (errno) {
               LOG_FATAL("%s\n", "errpr during int conversion");
               exit(1);
          }
          free(dp_char);
 
          if (dp_filter->min > 0 && cov < dp_filter->min) {
               vcf_var_add_to_filter(var, dp_filter->id_min);
          }
          if (dp_filter->max > 0 && cov > dp_filter->max) {
               vcf_var_add_to_filter(var, dp_filter->id_max);
          }
     }
}


void apply_snvqual_threshold(var_t *var, snvqual_filter_t *snvqual_filter)
{
     assert (! vcf_var_has_info_key(NULL, var, "INDEL"));
     if (! snvqual_filter->thresh) {
          return;
     }
     if (var->qual>-1 && var->qual<snvqual_filter->thresh) {
          vcf_var_add_to_filter(var, snvqual_filter->id);
     }
}


void apply_indelqual_threshold(var_t *var, indelqual_filter_t *indelqual_filter)
{
     assert (vcf_var_has_info_key(NULL, var, "INDEL"));
     if (! indelqual_filter->thresh) {
          return;
     }
     if (var->qual>-1 && var->qual<indelqual_filter->thresh) {
          vcf_var_add_to_filter(var, indelqual_filter->id);
     }
}


void apply_sb_threshold(var_t *var, sb_filter_t *sb_filter)
{
     char *sb_char = NULL;
     int sb;

     if (! sb_filter->thresh) {
          return;
     }

     if ( ! vcf_var_has_info_key(&sb_char, var, "SB")) {
          if ( ! sb_missing_warning_printed) {
               LOG_WARN("%s\n", "Requested SB filtering failed since SB tag is missing in variant");
               sb_missing_warning_printed = 1;
          }
          return;
     }
     sb = atoi(sb_char);
     free(sb_char);

     if (sb > sb_filter->thresh) {
          if (sb_filter->no_compound || alt_mostly_on_one_strand(var)) {
               vcf_var_add_to_filter(var, sb_filter->id);
          }
     }
}


/* returns -1 on error 
 *
 * filter everything that's not significant
 * 
 * Very similar to apply_sb_filter_mtc, but reverse testing logic and only looking at non consvars
 *
 * Will ignore indels
 */
int apply_snvqual_filter_mtc(snvqual_filter_t *snvqual_filter, var_t **vars, const long int num_vars)
{
     /* can only apply this logic to variants that are not consensus
      * variants, i.e those that actually have a quality. therefore
      * keep track of non cons var indeces */
     long int *orig_idx = NULL; /* of size num_noncons_vars */
     double *noncons_errprobs = NULL;
     long int num_noncons_vars = 0;
     long int i;

     if (snvqual_filter->ntests && num_vars > snvqual_filter->ntests) {
         LOG_WARN("%s\n", "Number of predefined tests for snvqual filter larger than number of variants! Are you sure that makes sense?");
     }

     /* collect values from noncons vars only and keep track of their indeces
      */
     orig_idx = malloc(num_vars * sizeof(long int));
     if ( ! orig_idx) {
          LOG_FATAL("%s\n", "out of memory");
          return -1;
     }
     noncons_errprobs = malloc(num_vars * sizeof(double));
     if ( ! noncons_errprobs) {
          LOG_FATAL("%s\n", "out of memory");
          return -1;
     }
     num_noncons_vars = 0;
     for (i=0; i<num_vars; i++) {
          if (vars[i]->qual>-1 && ! vcf_var_has_info_key(NULL, vars[i], "INDEL")) {
               noncons_errprobs[num_noncons_vars] = PHREDQUAL_TO_PROB(vars[i]->qual);
               orig_idx[num_noncons_vars] = i;
               num_noncons_vars += 1;
          }
     }
     if (! num_noncons_vars) {
          free(noncons_errprobs);
          free(orig_idx);
          return 0;
     }
     orig_idx = realloc(orig_idx, (num_noncons_vars * sizeof(long int)));
     noncons_errprobs = realloc(noncons_errprobs, (num_noncons_vars * sizeof(double)));

     /* only now we can set the number of tests (if it wasn't set by
      * caller) */
     if (! snvqual_filter->ntests) {
          snvqual_filter->ntests = num_noncons_vars;
     }

     /* multiple testing correction
      */
     if (snvqual_filter->mtc_type == MTC_BONF) {
          bonf_corr(noncons_errprobs, num_noncons_vars, 
                    snvqual_filter->ntests);
          
     } else if (snvqual_filter->mtc_type == MTC_HOLMBONF) {
          holm_bonf_corr(noncons_errprobs, num_noncons_vars, 
                         snvqual_filter->alpha, snvqual_filter->ntests);
          
     } else if (snvqual_filter->mtc_type == MTC_FDR) {
          long int num_rej = 0;
          long int *idx_rej; /* indices of rejected i.e. significant values */
          
          num_rej = fdr(noncons_errprobs, num_noncons_vars, 
                        snvqual_filter->alpha, snvqual_filter->ntests, 
                        &idx_rej);
          for (i=0; i<num_rej; i++) {
               long int idx = idx_rej[i];
               noncons_errprobs[idx] = -1;
          }
          free(idx_rej);
          
     } else {
          LOG_FATAL("Internal error: unknown MTC type %d\n", snvqual_filter->mtc_type);
          return -1;
     }
     
     for (i=0; i<num_noncons_vars; i++) {
          if (noncons_errprobs[i] > snvqual_filter->alpha) {
               vcf_var_add_to_filter(vars[orig_idx[i]], snvqual_filter->id);
          }
     }

     free(orig_idx);
     free(noncons_errprobs);

     return 0;
}



/* returns -1 on error 
 *
 * filter everything that's not significant
 * 
 * Very similar to apply_sb_filter_mtc, but reverse testing logic and only looking at non consvars
 *
 */
int apply_indelqual_filter_mtc(indelqual_filter_t *indelqual_filter, var_t **vars, const long int num_vars)
{
     /* can only apply this logic to variants that are not consensus
      * variants, i.e those that actually have a quality. therefore
      * keep track of non cons var indeces */
     long int *orig_idx = NULL; /* of size num_noncons_vars */
     double *noncons_errprobs = NULL;
     long int num_noncons_vars = 0;
     long int i;

     /* FIXME function almost identical to apply_indelqual_filter_mtc just different filter can be easily merged by accepting both types of variants */

     if (indelqual_filter->ntests && num_vars > indelqual_filter->ntests) {
         LOG_WARN("%s\n", "Number of predefined tests for indelqual filter larger than number of variants! Are you sure that makes sense?");
     }

     /* collect values from noncons vars only and keep track of their indeces
      */
     orig_idx = malloc(num_vars * sizeof(long int));
     if ( ! orig_idx) {
          LOG_FATAL("%s\n", "out of memory");
          return -1;
     }
     noncons_errprobs = malloc(num_vars * sizeof(double));
     if ( ! noncons_errprobs) {
          LOG_FATAL("%s\n", "out of memory");
          return -1;
     }
     num_noncons_vars = 0;
     for (i=0; i<num_vars; i++) {
          if (vars[i]->qual>-1 && vcf_var_has_info_key(NULL, vars[i], "INDEL")) {
               noncons_errprobs[num_noncons_vars] = PHREDQUAL_TO_PROB(vars[i]->qual);
               orig_idx[num_noncons_vars] = i;
               num_noncons_vars += 1;
          }
     }
     if (! num_noncons_vars) {
          free(noncons_errprobs);
          free(orig_idx);
          return 0;
     }
     orig_idx = realloc(orig_idx, (num_noncons_vars * sizeof(long int)));
     noncons_errprobs = realloc(noncons_errprobs, (num_noncons_vars * sizeof(double)));

     /* only now we can set the number of tests (if it wasn't set by
      * caller) */
     if (! indelqual_filter->ntests) {
          indelqual_filter->ntests = num_noncons_vars;
     }

     /* multiple testing correction
      */
     if (indelqual_filter->mtc_type == MTC_BONF) {
          bonf_corr(noncons_errprobs, num_noncons_vars, 
                    indelqual_filter->ntests);
          
     } else if (indelqual_filter->mtc_type == MTC_HOLMBONF) {
          holm_bonf_corr(noncons_errprobs, num_noncons_vars, 
                         indelqual_filter->alpha, indelqual_filter->ntests);
          
     } else if (indelqual_filter->mtc_type == MTC_FDR) {
          long int num_rej = 0;
          long int *idx_rej; /* indices of rejected i.e. significant values */
          
          num_rej = fdr(noncons_errprobs, num_noncons_vars, 
                        indelqual_filter->alpha, indelqual_filter->ntests, 
                        &idx_rej);
          for (i=0; i<num_rej; i++) {
               long int idx = idx_rej[i];
               noncons_errprobs[idx] = -1;
          }
          free(idx_rej);
          
     } else {
          LOG_FATAL("Internal error: unknown MTC type %d\n", indelqual_filter->mtc_type);
          free(orig_idx);
          free(noncons_errprobs);
          return -1;
     }
     
     for (i=0; i<num_noncons_vars; i++) {
          if (noncons_errprobs[i] > indelqual_filter->alpha) {
               vcf_var_add_to_filter(vars[orig_idx[i]], indelqual_filter->id);
          }
     }

     free(orig_idx);
     free(noncons_errprobs);

     return 0;
}


/* returns -1 on error 
 *
 * filter everything that's significant
 *
 * very similar to in apply_snvqual_filter_mtc, but reverse logic and looking at all vars
 */
int apply_sb_filter_mtc(sb_filter_t *sb_filter, var_t **vars, const long int num_vars)
{
     double *sb_probs = NULL;
     long int i;
     int num_ign = 0;
     long int *real_idx;/* we might ignore some variants (missing values etc). keep track of real indices of kept vars */

     
     /* collect values from vars kept in mem
      */
     sb_probs = malloc(num_vars * sizeof(double));
     if ( ! sb_probs) {LOG_FATAL("%s\n", "out of memory"); return -1;}
     real_idx = malloc(num_vars * sizeof(long int));
     if ( ! real_idx) {LOG_FATAL("%s\n", "out of memory"); return -1;}

     for (i=0; i<num_vars; i++) {
          char *sb_char = NULL;
          
          /* count indels if sb filter is not too be applied */
          if (! sb_filter->incl_indels && vcf_var_is_indel(vars[i])) {
               num_ign += 1;
               continue;
          }

          if ( ! vcf_var_has_info_key(&sb_char, vars[i], "SB")) {
               if ( ! sb_missing_warning_printed) {
                    LOG_WARN("%s\n", "At least one variant has no SB tag! SB filtering will be incomplete");
                    sb_missing_warning_printed = 1;
               }
               num_ign += 1;
               continue;
          }

          sb_probs[i-num_ign] = PHREDQUAL_TO_PROB(atoi(sb_char));
          real_idx[i-num_ign] = i;
          free(sb_char);
     }
     sb_probs = realloc(sb_probs, (num_vars-num_ign) * sizeof(double));
     real_idx = realloc(real_idx, (num_vars-num_ign) * sizeof(long int));

     if (! sb_filter->ntests) {
          sb_filter->ntests = num_vars - num_ign;
     } else {
          if (num_vars-num_ign > sb_filter->ntests) {
               LOG_WARN("%s\n", "Number of predefined tests for SB filter larger than number of variants! Are you sure that makes sense?");
          }
     }


     /* multiple testing correction
      */
     if (sb_filter->mtc_type == MTC_BONF) {
          bonf_corr(sb_probs, num_vars-num_ign, 
                    sb_filter->ntests);
          
     } else if (sb_filter->mtc_type == MTC_HOLMBONF) {
          holm_bonf_corr(sb_probs, num_vars-num_ign, 
                         sb_filter->alpha, sb_filter->ntests);
          
     } else if (sb_filter->mtc_type == MTC_FDR) {
          long int num_rej = 0;
          long int *idx_rej; /* indices of rejected i.e. significant values */
          
          num_rej = fdr(sb_probs, num_vars-num_ign, 
                        sb_filter->alpha, sb_filter->ntests, 
                        &idx_rej);
          for (i=0; i<num_rej; i++) {
               long int idx = idx_rej[i];
               sb_probs[real_idx[idx]] = -1;
          }
          free(idx_rej);
          
     } else {
          LOG_FATAL("Internal error: unknown MTC type %d\n", sb_filter->mtc_type);
          return -1;
     }
     
     for (i=0; i<num_vars; i++) {
          if (sb_probs[i] < sb_filter->alpha) {
               if (sb_filter->no_compound || alt_mostly_on_one_strand(vars[i])) {
                    vcf_var_add_to_filter(vars[i], sb_filter->id);
               }
          }
     }

     free(real_idx);
     free(sb_probs);

     return 0;
}


/* adds FILTER tags to vcf header based on config. also initializes
 * filter ids!
 */
void cfg_filter_to_vcf_header(filter_conf_t *cfg, char **header)
{
     char full_filter_str[FILTER_STRSIZE];

     /* for getting rid of all those trailing float zeros we might want to look at
        http://stackoverflow.com/questions/277772/avoid-trailing-zeroes-in-printf */

     if (cfg->af_filter.min > 0) {
          snprintf(cfg->af_filter.id_min, FILTER_ID_STRSIZE, "min_af_%f", cfg->af_filter.min);
          snprintf(full_filter_str, FILTER_STRSIZE,
                   "##FILTER=<ID=%s,Description=\"Minimum allele frequency %f\">\n",
                   cfg->af_filter.id_min, cfg->af_filter.min);
          vcf_header_add(header, full_filter_str);
     }

     if (cfg->af_filter.max > 0) {
          snprintf(cfg->af_filter.id_max, FILTER_ID_STRSIZE, "max_af_%f", cfg->af_filter.max);
          snprintf(full_filter_str, FILTER_STRSIZE,
               "##FILTER=<ID=%s,Description=\"Maximum allele frequency %f\">\n",
               cfg->af_filter.id_max, cfg->af_filter.max);
          vcf_header_add(header, full_filter_str);
     }

     if (cfg->dp_filter.min > 0) {
          snprintf(cfg->dp_filter.id_min, FILTER_ID_STRSIZE, "min_dp_%d", cfg->dp_filter.min);
          snprintf(full_filter_str, FILTER_STRSIZE,
               "##FILTER=<ID=%s,Description=\"Minimum Coverage %d\">\n",
               cfg->dp_filter.id_min, cfg->dp_filter.min);
          vcf_header_add(header, full_filter_str);
     }
     if (cfg->dp_filter.max > 0) {
          snprintf(cfg->dp_filter.id_max, FILTER_ID_STRSIZE, "max_dp_%d", cfg->dp_filter.max);
          snprintf(full_filter_str, FILTER_STRSIZE,
               "##FILTER=<ID=%s,Description=\"Maximum Coverage %d\">\n",
               cfg->dp_filter.id_max, cfg->dp_filter.max);
          vcf_header_add(header, full_filter_str);
     }

     assert (! (cfg->sb_filter.thresh > 0 && cfg->sb_filter.mtc_type != MTC_NONE));
     if (cfg->sb_filter.thresh > 0) {
          snprintf(cfg->sb_filter.id, FILTER_ID_STRSIZE, "max_sb_%d", cfg->sb_filter.thresh);
          snprintf(full_filter_str, FILTER_STRSIZE,
               "##FILTER=<ID=%s,Description=\"Maximum Strand-Bias Phred %d\">\n",
               cfg->sb_filter.id, cfg->sb_filter.thresh);
          vcf_header_add(header, full_filter_str);
          
     } else if (cfg->sb_filter.mtc_type != MTC_NONE) {
          char buf[64];
          mtc_str(buf, cfg->sb_filter.mtc_type);
          snprintf(cfg->sb_filter.id, FILTER_ID_STRSIZE, "sb_%s", buf);
          snprintf(full_filter_str, FILTER_STRSIZE,
               "##FILTER=<ID=%s,Description=\"Strand-Bias Multiple Testing Correction: %s corr. pvalue > %f\">\n",
                   cfg->sb_filter.id, buf, cfg->sb_filter.alpha);
          vcf_header_add(header, full_filter_str);
     }

     assert (! (cfg->snvqual_filter.thresh > 0 && cfg->snvqual_filter.mtc_type != MTC_NONE));
     if (cfg->snvqual_filter.thresh > 0) {
          snprintf(cfg->snvqual_filter.id, FILTER_ID_STRSIZE, "min_snvqual_%d", cfg->snvqual_filter.thresh);
          snprintf(full_filter_str, FILTER_STRSIZE,
               "##FILTER=<ID=%s,Description=\"Minimum SNV Quality (Phred) %d\">\n",
               cfg->snvqual_filter.id, cfg->snvqual_filter.thresh);
          vcf_header_add(header, full_filter_str);
          
     } else if (cfg->snvqual_filter.mtc_type != MTC_NONE) {
          char buf[64];
          mtc_str(buf, cfg->snvqual_filter.mtc_type);
          snprintf(cfg->snvqual_filter.id, FILTER_ID_STRSIZE, "snvqual_%s", buf);
          snprintf(full_filter_str, FILTER_STRSIZE,
               "##FILTER=<ID=%s,Description=\"SNV Quality Multiple Testing Correction: %s corr. pvalue < %f\">\n",
                   cfg->snvqual_filter.id, buf, cfg->snvqual_filter.alpha);
          vcf_header_add(header, full_filter_str);
     }

     assert (! (cfg->indelqual_filter.thresh > 0 && cfg->indelqual_filter.mtc_type != MTC_NONE));
     if (cfg->indelqual_filter.thresh > 0) {
          snprintf(cfg->indelqual_filter.id, FILTER_ID_STRSIZE, "min_indelqual_%d", cfg->indelqual_filter.thresh);
          snprintf(full_filter_str, FILTER_STRSIZE,
               "##FILTER=<ID=%s,Description=\"Minimum Indel Quality (Phred) %d\">\n",
               cfg->indelqual_filter.id, cfg->indelqual_filter.thresh);
          vcf_header_add(header, full_filter_str);
          
     } else if (cfg->indelqual_filter.mtc_type != MTC_NONE) {
          char buf[64];
          mtc_str(buf, cfg->indelqual_filter.mtc_type);
          snprintf(cfg->indelqual_filter.id, FILTER_ID_STRSIZE, "indelqual_%s", buf);
          snprintf(full_filter_str, FILTER_STRSIZE,
               "##FILTER=<ID=%s,Description=\"Indel Quality Multiple Testing Correction: %s corr. pvalue < %f\">\n",
                   cfg->indelqual_filter.id, buf, cfg->indelqual_filter.alpha);
          vcf_header_add(header, full_filter_str);
     }
}


int
main_filter(int argc, char *argv[])
{
     filter_conf_t cfg;
     char *vcf_in = NULL, *vcf_out = NULL;
     static int print_only_passed = 1;
     static int sb_filter_no_compound = 0;
     static int sb_filter_incl_indels = 0;
     static int only_indels = 0;
     static int only_snvs = 0;
     char *vcf_header = NULL;
     var_t **vars = NULL;
     long int num_vars = 0; /* isn't long overkill here ? */
     long int vars_size = 0; /* keeping track of how much memory we've got pre-allocated */
     long int i;
     static int no_defaults = 0;

     /* default filter options */
     memset(&cfg, 0, sizeof(filter_conf_t));
     cfg.dp_filter.min = cfg.dp_filter.max = -1;
     cfg.af_filter.min = cfg.af_filter.max = -1;
     cfg.sb_filter.alpha = DEFAULT_SIG;
     cfg.snvqual_filter.alpha = DEFAULT_SIG;
     cfg.indelqual_filter.alpha = DEFAULT_SIG;


    /* keep in sync with long_opts_str and usage
     *
     * getopt is a pain in the whole when it comes to syncing of long
     * and short args and usage. check out gopt, libcfu...
     */
    while (1) {
         int c;
         static struct option long_opts[] = {
              /* see usage sync */
              {"verbose", no_argument, &verbose, 1},
              {"debug", no_argument, &debug, 1},
              {"print-all", no_argument, &print_only_passed, 0},
              {"no-defaults", no_argument, &no_defaults, 1},
              {"only-indels", no_argument, &only_indels, 1},
              {"only-snvs", no_argument, &only_snvs, 1},

              {"help", no_argument, NULL, 'h'},
              {"in", required_argument, NULL, 'i'},
              {"out", required_argument, NULL, 'o'},

              {"cov-min", required_argument, NULL, 'v'},
              {"cov-max", required_argument, NULL, 'V'},

              {"af-min", required_argument, NULL, 'a'},
              {"af-max", required_argument, NULL, 'A'},

              {"sb-thresh", required_argument, NULL, 'B'},
              {"sb-mtc", required_argument, NULL, 'b'},
              {"sb-alpha", required_argument, NULL, 'c'},
              {"sb-no-compound", no_argument, &sb_filter_no_compound, 1},
              {"sb-incl-indels", no_argument, &sb_filter_incl_indels, 1},

              {"snvqual-thresh", required_argument, NULL, 'Q'},
              {"snvqual-mtc", required_argument, NULL, 'q'},
              {"snvqual-alpha", required_argument, NULL, 'r'},
              {"snvqual-ntests", required_argument, NULL, 's'},

              {"indelqual-thresh", required_argument, NULL, 'K'},
              {"indelqual-mtc", required_argument, NULL, 'k'},
              {"indelqual-alpha", required_argument, NULL, 'l'},
              {"indelqual-ntests", required_argument, NULL, 'm'},

              {0, 0, 0, 0} /* sentinel */
         };

         /* keep in sync with long_opts and usage */
         static const char *long_opts_str = "hi:o:v:V:a:A:B:b:c:Q:q:r:s:K:k:l:m:";

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
              usage(& cfg);
              return 0;

         case 'i':
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

         case 'v':
              if (! isdigit(optarg[0])) {
                   LOG_FATAL("Non-numeric argument provided: %s\n", optarg);
                   return -1;
              }
              cfg.dp_filter.min = atoi(optarg);
              break;
         case 'V':
              if (! isdigit(optarg[0])) {
                   LOG_FATAL("Non-numeric argument provided: %s\n", optarg);
                   return -1;
              }
              cfg.dp_filter.max = atoi(optarg);
              break;

         case 'a':
              if (! isdigit(optarg[0])) {
                   LOG_FATAL("Non-numeric argument provided: %s\n", optarg);
                   return -1;
              }
              cfg.af_filter.min = strtof(optarg, NULL);
              break;
         case 'A':
              if (! isdigit(optarg[0])) {
                   LOG_FATAL("Non-numeric argument provided: %s\n", optarg);
                   return -1;
              }
              cfg.af_filter.max = strtof(optarg, NULL);
              break;

         case 'B':
              if (! isdigit(optarg[0])) {
                   LOG_FATAL("Non-numeric argument provided: %s\n", optarg);
                   return -1;
              }
              cfg.sb_filter.thresh = atoi(optarg);
              break;
         case 'b':
              cfg.sb_filter.mtc_type = mtc_str_to_type(optarg);
              if (-1 == cfg.sb_filter.mtc_type) {
                   LOG_FATAL("Unknown multiple testing correction type '%s' for strandbias filtering\n", optarg);
                   return -1;
              }
              break;
         case 'c':
              if (! isdigit(optarg[0])) {
                   LOG_FATAL("Non-numeric argument provided: %s\n", optarg);
                   return -1;
              }
              cfg.sb_filter.alpha = strtof(optarg, NULL);
              break;

         case 'Q':
              if (! isdigit(optarg[0])) {
                   LOG_FATAL("Non-numeric argument provided: %s\n", optarg);
                   return -1;
              }
              cfg.snvqual_filter.thresh = atoi(optarg);
              break;
         case 'q':
              cfg.snvqual_filter.mtc_type = mtc_str_to_type(optarg);
              if (-1 == cfg.snvqual_filter.mtc_type) {
                   LOG_FATAL("Unknown multiple testing correction type '%s' for snv quality filtering\n", optarg);
                   return -1;
              }
              break;
         case 'r':
              if (! isdigit(optarg[0])) {
                   LOG_FATAL("Non-numeric argument provided: %s\n", optarg);
                   return -1;
              }
              cfg.snvqual_filter.alpha = strtof(optarg, NULL);
              break;
         case 's':
              if (! isdigit(optarg[0])) {
                   LOG_FATAL("Non-numeric argument provided: %s\n", optarg);
                   return -1;
              }
              cfg.snvqual_filter.ntests = atol(optarg);
              break;

         case 'K':
              if (! isdigit(optarg[0])) {
                   LOG_FATAL("Non-numeric argument provided: %s\n", optarg);
                   return -1;
              }
              cfg.indelqual_filter.thresh = atoi(optarg);
              break;
         case 'k':
              cfg.indelqual_filter.mtc_type = mtc_str_to_type(optarg);
              if (-1 == cfg.indelqual_filter.mtc_type) {
                   LOG_FATAL("Unknown multiple testing correction type '%s' for snv quality filtering\n", optarg);
                   return -1;
              }
              break;
         case 'l':
              if (! isdigit(optarg[0])) {
                   LOG_FATAL("Non-numeric argument provided: %s\n", optarg);
                   return -1;
              }
              cfg.indelqual_filter.alpha = strtof(optarg, NULL);
              break;
         case 'm':
              if (! isdigit(optarg[0])) {
                   LOG_FATAL("Non-numeric argument provided: %s\n", optarg);
                   return -1;
              }
              cfg.indelqual_filter.ntests = atol(optarg);
              break;

         case '?':
              LOG_FATAL("%s\n", "Unrecognized argument found. Exiting...\n");
              return 1;

         default:
              break;
         }
    }
    cfg.print_only_passed = print_only_passed;
    cfg.only_indels = only_indels;
    cfg.only_snvs = only_snvs;
    cfg.sb_filter.no_compound = sb_filter_no_compound;
    cfg.sb_filter.incl_indels = sb_filter_incl_indels;

    if (cfg.only_indels && cfg.only_snvs) {
         LOG_FATAL("%s\n", "Can't keep only indels and only snvs");
         return 1;
    }
    
    if (! no_defaults) {
         if (cfg.sb_filter.mtc_type==MTC_NONE && ! cfg.sb_filter.thresh) {
              LOG_VERBOSE("%s\n", "Setting default SB filtering method to holm-bonf");
              cfg.sb_filter.mtc_type = MTC_FDR;
              cfg.sb_filter.alpha = 0.001;
         }
         if (cfg.dp_filter.min<0) {
              LOG_VERBOSE("%s\n", "Setting default minimum coverage to 10");
              cfg.dp_filter.min = 10;
         }
    } else {
         LOG_VERBOSE("%s\n", "Skipping default settings");
    }

    if (0 != argc - optind - 1) {/* FIXME needed at all? */
         LOG_FATAL("%s\n", "Unrecognized argument found. Exiting...\n");
         return 1;
    }

    /* logic check of command line parameters
     */
    if (cfg.dp_filter.max > 0 &&  cfg.dp_filter.max < cfg.dp_filter.min) {
         LOG_FATAL("%s\n", "Invalid coverage-filter settings");
         return 1;
    }
    if ((cfg.af_filter.max > 0 && cfg.af_filter.max < cfg.af_filter.min) ||
        (cfg.af_filter.max > 1.0)) {
         LOG_FATAL("%s\n", "Invalid AF-filter settings");
         return 1;
    }

    if (cfg.sb_filter.thresh && cfg.sb_filter.mtc_type != MTC_NONE) {
         LOG_FATAL("%s\n", "Can't use fixed strand-bias threshold *and* multiple testing correction.");
         return 1;
    }
    if (cfg.snvqual_filter.thresh && cfg.snvqual_filter.mtc_type != MTC_NONE) {
         LOG_FATAL("%s\n", "Can't use fixed SNV quality threshold *and* multiple testing correction.");
         return 1;
    }
    if (cfg.indelqual_filter.thresh && cfg.indelqual_filter.mtc_type != MTC_NONE) {
         LOG_FATAL("%s\n", "Can't use fixed indel quality threshold *and* multiple testing correction.");
         return 1;
    }

    if (argc == 2) {
        fprintf(stderr, "\n");
        usage(& cfg);
        return 1;
    }

    if (debug) {
          dump_filter_conf(& cfg);
     }

    /* missing file args default to stdin and stdout
     */
    if  (! vcf_in) {
         vcf_in = malloc(2 * sizeof(char));
         strcpy(vcf_in, "-");
    }
    if  (! vcf_out) {
         vcf_out = malloc(2 * sizeof(char));
         strcpy(vcf_out, "-");
    }
    LOG_DEBUG("vcf_in=%s vcf_out=%s\n", vcf_in, vcf_out);


    /* open vcf files
     */
    if (vcf_file_open(& cfg.vcf_in, vcf_in,
                      HAS_GZIP_EXT(vcf_in), 'r')) {
         LOG_ERROR("Couldn't open %s\n", vcf_in);
         return 1;
    }
    if (vcf_file_open(& cfg.vcf_out, vcf_out,
                      HAS_GZIP_EXT(vcf_out), 'w')) {
         LOG_ERROR("Couldn't open %s\n", vcf_out);
         return 1;
    }
    free(vcf_in);
    free(vcf_out);

    /* FIXME everything below here should go into a function with args:
       - cfg
       - ...what else?
    */

    /* print header
     */
    if (0 !=  vcf_parse_header(&vcf_header, & cfg.vcf_in)) {
         /* LOG_WARN("%s\n", "vcf_parse_header() failed"); */
         if (vcf_file_seek(& cfg.vcf_in, 0, SEEK_SET)) {
              LOG_FATAL("%s\n", "Couldn't rewind file to parse variants"
                        " after header parsing failed");
              return -1;
         }
    }
    /* also sets filter names */
    cfg_filter_to_vcf_header(& cfg, &vcf_header);
    vcf_write_header(& cfg.vcf_out, vcf_header);
    free(vcf_header);


    /* read in variants. since many filters perform multiple testing
     * correction and therefore need to look at all variants we keep
     * it simple and load them all into memory. 
     * 
     * in theory we could apply all 'simple' filters directly within
     * the loop here and depending on the result spit the variant out
     * or not. only complex filters need to see all variants first to,
     * e.g. apply multiple testing.
     */
    num_vars = 0;
    while (1) {
         var_t *var;
         int rc;
         int is_indel = 0;

         vcf_new_var(&var);
         rc = vcf_parse_var(& cfg.vcf_in, var);
         if (rc) {
              /* how to distinguish between error and EOF? */
              free(var);
              break;
         }

         is_indel = vcf_var_is_indel(var);

         if (cfg.only_snvs && is_indel) {
              free(var);
              continue;
         } else if (cfg.only_indels && ! is_indel) {
              free(var);
              continue;
         }

         /* read all in, no matter if already filtered. we keep adding filters */
         num_vars +=1;
         if (num_vars >= vars_size) {
              const long incr = 128;
              vars = realloc(vars, (vars_size+incr) * sizeof(var_t*));
              vars_size += incr;
         }
         vars[num_vars-1] = var;
#ifdef TRACE
         {
              char *key;
              vcf_var_key(&key,  vars[num_vars-1]);
              fprintf(stderr, "storing var %ld+1: %s\n", num_vars, key);
              free(key);
         }
#endif

         /* filters applying to all types of variants
          */
         apply_af_filter(var, & cfg.af_filter);
         apply_dp_filter(var, & cfg.dp_filter);

         /* quality threshold per variant type
          */
         if (! is_indel) {
              if (cfg.snvqual_filter.thresh) {
                   assert(cfg.snvqual_filter.mtc_type == MTC_NONE);
                   apply_snvqual_threshold(var, & cfg.snvqual_filter);
              }

         } else {
              if (cfg.indelqual_filter.thresh) {
                   assert(cfg.indelqual_filter.mtc_type == MTC_NONE);
                   apply_indelqual_threshold(var, & cfg.indelqual_filter);
              }
         }
         
         if (cfg.sb_filter.thresh) {
              if (! is_indel || cfg.sb_filter.incl_indels) {
                   assert(cfg.sb_filter.mtc_type == MTC_NONE);
                   apply_sb_threshold(var, & cfg.sb_filter);
              }
         }
    }

    if (num_vars) {
         vars = realloc(vars, (num_vars * sizeof(var_t*)));
    }
    vcf_file_close(& cfg.vcf_in);
    LOG_VERBOSE("Parsed %d variants\n", num_vars);


    if (cfg.sb_filter.mtc_type != MTC_NONE) {
         if (apply_sb_filter_mtc(& cfg.sb_filter, vars, num_vars)) {
              LOG_FATAL("%s\n", "Multiple testing correction on strand-bias pvalues failed");
              return -1;
         }
    }

    if (cfg.snvqual_filter.mtc_type != MTC_NONE) {
         if (apply_snvqual_filter_mtc(& cfg.snvqual_filter, vars, num_vars)) {
              LOG_FATAL("%s\n", "Multiple testing correction on SNV qualities failed");
              return -1;
         }
    }

    if (cfg.indelqual_filter.mtc_type != MTC_NONE) {
         if (apply_indelqual_filter_mtc(& cfg.indelqual_filter, vars, num_vars)) {
              LOG_FATAL("%s\n", "Multiple testing correction on Indel qualities failed");
              return -1;
         }
    }

    /* output
     */
    for (i=0; i<num_vars; i++) {
         var_t *v = vars[i];

         if (cfg.print_only_passed && ! (VCF_VAR_PASSES(v))) {
              continue;
         }

         /* add pass if no filters were set */
         if (! v->filter || strlen(v->filter)<=1) {
              char pass_str[] = "PASS";
              if (v->filter) {
                   free(v->filter);
              }
              v->filter = strdup(pass_str);
         }

         vcf_write_var(& cfg.vcf_out, v);
    }
    vcf_file_close(& cfg.vcf_out);


    for (i=0; i<num_vars; i++) {
         vcf_free_var(& vars[i]);
    }
    free(vars);

    LOG_VERBOSE("%s\n", "Successful exit.");

    return 0;
}
/* main_filter */


/* gcc lofreq_filter.c -o lofreq_filter -I../lofreq_core -I../uthash/ ../lofreq_core/liblofreq_core.a   -lz -DMAIN_FILTER */
#ifdef MAIN_FILTER

int
main(int argc, char *argv[])
{
     return main_filter(argc+1, argv-1);
}
#endif
