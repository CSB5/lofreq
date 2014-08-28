/* -*- c-file-style: "k&r"; indent-tabs-mode: nil; -*- */


/*********************************************************************
 *
 * FIXME update license
 *
 *********************************************************************/

/*
 * This file is partially based on samtools' bam_plcmd.c Whenever
 * parts of the code looks like they've been written by a other-wordly
 * wizard, then it was probably from Heng Li.
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

/* libbam includes */
#include "faidx.h"
#include "sam.h"
#include "kstring.h"
/* from bedidx.c */
void *bed_read(const char *fn);
void bed_destroy(void *_h);
int bed_overlap(const void *_h, const char *chr, int beg, int end);

/* lofreq includes */
#include "snpcaller.h"
#include "vcf.h"
#include "fet.h"
#include "utils.h"
#include "log.h"
#include "plp.h"
#include "bed.h"
#include "defaults.h"

#if 1
#define MYNAME "lofreq call"
#else
#define MYNAME PACKAGE
#endif


#define BUF_SIZE 1<<16



/* number of tests performed (CONSVAR doesn't count). for downstream
 * multiple testing correction. corresponds to bonf if bonf_dynamic is
 * true. */
long long int num_snv_tests = 0;
long long int num_indel_tests = 0;
/* FIXME extend to keep some more stats, e.g. num_pos_with_cov etc */




void
report_var(vcf_file_t *vcf_file, const plp_col_t *p, const char *ref, 
           const char *alt, const float af, const int qual,
           const int is_indel, const int is_consvar, 
           const dp4_counts_t *dp4)
{
     var_t *var;
     double sb_left_pv, sb_right_pv, sb_two_pv;
     int sb_qual;
#if MQ_BIAS_FIXME
     int *ref_mq = NULL, *alt_mq = NULL;
     long int num_ref;
     long int num_alt;
     int ref_nt4 = bam_nt4_table[(int)ref];
     int alt_nt4 = bam_nt4_table[(int)alt];
#endif

     vcf_new_var(&var);
     var->chrom = strdup(p->target);
     var->pos = p->pos;
     /* var->id = NA */
     var->ref = strdup(ref);
     var->alt = strdup(alt);
     if (qual>-1) {
          var->qual = qual;
     }
     /* var->filter = NA */ 

#if MQ_BIAS_FIXME
     /* FIXME: this will not work for indel calling...what are these being used for? */
     assert (p->map_quals[ref_nt4].n == p->fw_counts[ref_nt4] + p->rv_counts[ref_nt4]);
     assert (p->map_quals[alt_nt4].n == p->fw_counts[alt_nt4] + p->rv_counts[alt_nt4]);
     num_alt = p->map_quals[alt_nt4].n;
     num_ref = p->map_quals[ref_nt4].n;
#endif
#if 0
     LOG_FIXME("%s\n", "fisher's test on unfiltered counts. Is that wanted?");
     LOG_FIXME("%s\n", "Hand down filtered quals/counts might solve fisher and mq bias problems. DP4 is supposed to be filtered version anyway");
#endif
     /* strand bias
      */
     /* double sb_prob = kt... Assignment removed to shut up clang static analyzer */
     (void) kt_fisher_exact(dp4->ref_fw, dp4->ref_rv, 
                            dp4->alt_fw, dp4->alt_rv,
                            &sb_left_pv, &sb_right_pv, &sb_two_pv);
     sb_qual = PROB_TO_PHREDQUAL_SAFE(sb_two_pv);

#if MQ_BIAS_FIXME
     /* mq bias
      */
     assert (p->map_quals[ref_nt4].n == p->base_quals[ref_nt4].n);
     assert (p->map_quals[alt_nt4].n == p->base_quals[alt_nt4].n);
     /* FIXME intro of n_alt and n_ref would make things easier here */
     if (NULL == (ref_mq = malloc(num_ref * sizeof(int)))) {
          /* coverage = base-count after read level filtering */
          fprintf(stderr, "FATAL: couldn't allocate memory at %s:%s():%d\n",
                  __FILE__, __FUNCTION__, __LINE__);
          return;
     }
     if (NULL == (alt_mq = malloc(num_alt * sizeof(int)))) {
          /* coverage = base-count after read level filtering */
          fprintf(stderr, "FATAL: couldn't allocate memory at %s:%s():%d\n",
                  __FILE__, __FUNCTION__, __LINE__);
          return;
     }
/*
     for (i=0; i<p->map_quals[ref_nt4].n; j++) {
          mq = p->map_quals[i].data[j];
     }
     for (i=0; i<p->map_quals[alt_nt4].n; j++) {
          mq = p->map_quals[i].data[j];
     }
*/
/* FIXME mq bias
   see mann_whitney_1947 in samtools
   need NUM_NT4 altbase
   need NUM_NT4 refbase

   LOG_FIXME("%s\n", "mq and bq filtering should be done here but needs conf");
   
   for (i=0; i<NUM_NT4; i++) {
     int nt = bam_nt4_rev_table[i];
     for (j=0; j<p->map_quals[i].n; j++) {
          mq = p->map_quals[i].data[j];
     }
}
*/
#endif
     int cover = p->coverage - p->num_tails;
	vcf_var_sprintf_info(var, is_indel? cover : p->coverage, af, sb_qual,
                          dp4, is_indel, is_consvar);

     vcf_write_var(vcf_file, var);
     vcf_free_var(&var);

#if MQ_BIAS_FIXME
     free(ref_mq);
     free(alt_mq);
#endif
}
/* report_var() */


/* report consensus substitution */
void
report_cons_sub(const plp_col_t *p, snvcall_conf_t *conf){
     
     const int is_indel = 0;
     const int is_consvar = 1;
     const int qual = -1;
     float af = base_count(p, p->cons_base[0]) / (float)p->coverage;

     char report_ref[2];
     report_ref[0] = p->ref_base;
     report_ref[1] = '\0';

     int ref_nt4;
     int alt_nt4;
     ref_nt4 = bam_nt4_table[(int)report_ref[0]];
     alt_nt4 = bam_nt4_table[(int)p->cons_base[0]];

     dp4_counts_t dp4;
     dp4.ref_fw = p->fw_counts[ref_nt4];
     dp4.ref_rv = p->rv_counts[ref_nt4];
     dp4.alt_fw = p->fw_counts[alt_nt4];
     dp4.alt_rv = p->rv_counts[alt_nt4];

     LOG_DEBUG("cons var snp: %s %d %c>%c\n",
               p->target, p->pos+1, p->ref_base, p->cons_base);          
     report_var(& conf->vcf_out, p, report_ref, p->cons_base, 
                af, qual, is_indel, is_consvar, &dp4);
}

/* report consensus insertion */
void
report_cons_ins(const plp_col_t *p, snvcall_conf_t *conf) {

     const int is_indel = 1;
     const int is_consvar = 1;
     const int qual = -1;

     char cons_ins_key[MAX_INDELSIZE];
     strcpy(cons_ins_key, p->cons_base+1);
     ins_event *it_ins = find_ins_sequence(&p->ins_event_counts, cons_ins_key);

     char report_ins_ref[2];
     char report_ins_alt[MAX_INDELSIZE];
     int ins_length = strlen(cons_ins_key);
     report_ins_ref[0] = report_ins_alt[0] = p->ref_base;
     int j;
     for (j = 0; j <= ins_length; ++j) {
          report_ins_alt[j+1] = cons_ins_key[j];
     }
     report_ins_ref[1] = report_ins_alt[j+1] = '\0';
     
     float af = it_ins->count / (float)p->coverage_indel;

     dp4_counts_t dp4;
     dp4.ref_fw = p->non_ins_fw_rv[0];
     dp4.ref_rv = p->non_ins_fw_rv[1];
     dp4.alt_fw = it_ins->fw_rv[0];
     dp4.alt_rv = it_ins->fw_rv[1];

     LOG_DEBUG("Consensus insertion: %s %d %s>%s\n",
               p->target, p->pos+1, report_ins_ref, report_ins_alt);
     report_var(& conf->vcf_out, p, report_ins_ref, report_ins_alt, 
                af, qual, is_indel, is_consvar, &dp4);
     return;
}

/* report consensus deletion */
void 
report_cons_del(const plp_col_t *p, snvcall_conf_t *conf) {
     
     const int is_indel = 1;
     const int is_consvar = 1;
     const int qual = -1;

     char cons_del_key[MAX_INDELSIZE];
     strcpy(cons_del_key, p->cons_base+1);
     del_event *it_del = find_del_sequence(&p->del_event_counts, cons_del_key);

     char report_del_ref[MAX_INDELSIZE];
     char report_del_alt[2];
     int del_length = strlen(cons_del_key);
     report_del_ref[0] = report_del_alt[0] = p->ref_base;
     int j;
     for (j = 0; j <= del_length; ++j) {
          report_del_ref[j+1] = cons_del_key[j];
     }
     report_del_ref[j+1] = report_del_alt[1] = '\0';

     float af = it_del->count / (float)p->coverage_indel;
     
     dp4_counts_t dp4;
     dp4.ref_fw = p->non_del_fw_rv[0];
     dp4.ref_rv = p->non_del_fw_rv[1];
     dp4.alt_fw = it_del->fw_rv[0];
     dp4.alt_rv = it_del->fw_rv[1];

     LOG_DEBUG("Consensus deletion: %s %d %s>%s\n",
               p->target, p->pos+1, report_del_ref, report_del_alt);
     report_var(&conf->vcf_out, p, report_del_ref, report_del_alt,
                af, qual, is_indel, is_consvar, &dp4);

}

int
call_alt_ins(const plp_col_t *p, double *bi_err_probs, int bi_num_err_probs,
             snvcall_conf_t *conf, ins_event *it) {

     int ins_counts[3];
     long double bi_pvalues[3];

     // prep for snpcaller
     ins_counts[0] = it->count;
     ins_counts[1] = ins_counts[2] = 0;
     LOG_DEBUG("%s %d: passing down %d quals with noncons_ins_counts"
               "(%d, %d, %d) to snpcaller()\n", p->target, p->pos+1,
               bi_num_err_probs, ins_counts[0], ins_counts[1], ins_counts[2]);
     // compute p-value for insertion
     if (snpcaller(bi_pvalues, bi_err_probs, bi_num_err_probs, ins_counts,
                   conf->bonf_indel, conf->sig)) {
          fprintf(stderr, "FATAL: snpcaller() failed at %s:%s():%d\n",
                  __FILE__, __FUNCTION__, __LINE__);
          return 1;
     }
     // see if there was an insertion
     long double bi_pvalue = bi_pvalues[0];
     if (bi_pvalue*conf->bonf_indel < conf->sig) {
          const int is_indel = 1;
          const int is_consvar = 0;
          const int qual = PROB_TO_PHREDQUAL(bi_pvalue);

          char report_ins_ref[2];
          char report_ins_alt[256];
          int ins_length = strlen(it->key);
          report_ins_ref[0] = report_ins_alt[0] = p->ref_base;
          int j;
          for (j = 0; j <= ins_length; ++j) {
               report_ins_alt[j+1] = it->key[j];
          }
          report_ins_ref[1] = report_ins_alt[j+1] = '\0';
          float af = it->count / (float)p->coverage_indel;
     
          dp4_counts_t dp4;
          dp4.ref_fw = p->non_ins_fw_rv[0];
          dp4.ref_rv = p->non_ins_fw_rv[1];
          dp4.alt_fw = it->fw_rv[0];
          dp4.alt_rv = it->fw_rv[1];

          LOG_DEBUG("Low freq insertion: %s %d %s>%s pv-prob:%g;pv-qual:%d\n",
                    p->target, p->pos+1, report_ins_ref, report_ins_alt,
                    bi_pvalue, qual);
          report_var(&conf->vcf_out, p, report_ins_ref, report_ins_alt,
                     af, qual, is_indel, is_consvar, &dp4);
     }
     return 0;
}

int call_alt_del(const plp_col_t *p, double *bd_err_probs, int bd_num_err_probs,
                 snvcall_conf_t *conf, del_event *it) {

     int del_counts[3];
     long double bd_pvalues[3];

     // prep for snpcaller
     del_counts[0] = it->count;
     del_counts[1] = del_counts[2] = 0;
     //int k;
     //for (k = 0; k < bd_num_err_probs; k++) {
     //     LOG_DEBUG("bd_err_prob: %lg\n", bd_err_probs[k]);
     //}
     LOG_DEBUG("%s %d: passing down %d quals with noncons_del_counts"
               "(%d, %d, %d) to snpcaller()\n", p->target, p->pos+1,
               bd_num_err_probs, del_counts[0], del_counts[1], del_counts[2]);

     // snpcaller for deletion
     if (snpcaller(bd_pvalues, bd_err_probs, bd_num_err_probs, del_counts,
                   conf->bonf_indel, conf->sig)) {
          fprintf(stderr, "FATAL: snpcaller() failed at %s:%s():%d\n",
                  __FILE__, __FUNCTION__, __LINE__);
          return 1;
     }
     // compute p-value deletion
     long double bd_pvalue = bd_pvalues[0];
     if (bd_pvalue*conf->bonf_indel < conf->sig) {
          const int is_indel = 1;
          const int is_consvar = 0;
          const int qual = PROB_TO_PHREDQUAL(bd_pvalue);

          char report_del_ref[256];
          char report_del_alt[2];
          int del_length = strlen(it->key);
          report_del_ref[0] = report_del_alt[0] = p->ref_base;
          int j;
          for (j = 0; j <= del_length; ++j) {
               report_del_ref[j+1] = it->key[j];
          }
          report_del_ref[j+1] = report_del_alt[1] = '\0';

          float af = it->count / (float)p->coverage_indel;
          
          dp4_counts_t dp4;
          dp4.ref_fw = p->non_del_fw_rv[0];
          dp4.ref_rv = p->non_del_fw_rv[1];
          dp4.alt_fw = it->fw_rv[0];
          dp4.alt_rv = it->fw_rv[1];

          LOG_DEBUG("Low freq deletion: %s %d %s>%s pv-prob:%g;pv-qual:%d\n",
                    p->target, p->pos+1, report_del_ref, report_del_alt,
                    bd_pvalue, qual);
          report_var(&conf->vcf_out, p, report_del_ref, report_del_alt,
                     af, qual, is_indel, is_consvar, &dp4);
     }
     return 0;
}

/* allocates bc_err_probs (to size bc_num_err_probs; also set here) and sets
 * values. user must free.
 *
 * qualities are merged here and filtering also happens here
 *
 * alt_bases, alt_counts and alt_raw_counts must be pre-allocated and
 * of size 3 and values will be set here (FIXME that makes the
 * function call awkward)
 */
void
plp_summary(const plp_col_t *plp_col, void* confp) 
{
     FILE* stream = stdout;
     snvcall_conf_t *conf = (snvcall_conf_t *)confp;
     static const char* title[] = {"BQ", "BAQ", "MQ", "SQ"};
     int i, x;

     fprintf(stream, "%s\t%d\t%c\t%s", plp_col->target, plp_col->pos+1,
             plp_col->ref_base, plp_col->cons_base);
     for (i=0; i<NUM_NT4; i++) {
          fprintf(stream, "\t%c:%lu/%lu",
                  bam_nt4_rev_table[i],
                  plp_col->fw_counts[i],
                  plp_col->rv_counts[i]);
     }

     fprintf(stream, "\theads:%d\ttails:%d", plp_col->num_heads, 
             plp_col->num_tails);
     fprintf(stream, "\tins=%d\tdels=%d", plp_col->num_ins, 
             plp_col->num_dels);
     fprintf(stream, "\n");
     for (i=0; i<NUM_NT4; i++) {
          int X = 3;/* bq, baq, mq */
#ifdef USE_SOURCEQUAL                         
          if (conf->flag & SNVCALL_USE_SQ) {
               X = 4;/* bq, baq, mq, sq */
          }
#endif
          for (x=0; x<X; x++) {
               int j;
               int nt = bam_nt4_rev_table[i];
               fprintf(stream, "  %s %c =", title[x], nt);
               for (j=0; j<plp_col->base_quals[i].n; j++) {
                    int q = -1;
                    if (x==0) {
                         q = plp_col->base_quals[i].data[j];
                    } else if (x==1 && conf->flag & SNVCALL_USE_BAQ) {
                         q = plp_col->baq_quals[i].data[j];
                    } else if (x==2) {
                         q = plp_col->map_quals[i].data[j];
                    } else if (x==3) {
                         q = plp_col->source_quals[i].data[j];
                    }
                    fprintf(stream, " %d", q);
               }
               fprintf(stream, "\n");
          }
     }
}


void
warn_old_fai(const char *fa)
{
     char *fai;
     if (!fa || fa[0]=='\0') {
          return;
     }

     fai = (char*) calloc(strlen(fa) + 5, 1);
     sprintf(fai, "%s.fai", fa);
     if (is_newer(fa, fai)==1) {
          LOG_WARN("Index for fasta file (%s) is older than fasta file! You should reindex (using faidx)!\n", fai);
     }
     free(fai);
}


/* low-freq vars always called against cons_base, which might be
 * different from ref_base. if cons_base != ref_base then it's a
 * cons-var.
 * 
 * Assuming conf->min_bq and read-level filtering was already done
 * upstream. altbase mangling happens here however.
 * 
 * FIXME split func into snv and indel calling
 */
void
call_vars(const plp_col_t *p, void *confp)
{
     double *bc_err_probs; /* error probs (qualities) passed down to snpcaller */
     int bc_num_err_probs; /* #elements in bc_err_probs */
     double *bi_err_probs, *bd_err_probs; /* error probs for indel calling */
     int bi_num_err_probs, bd_num_err_probs;
     int i;

     snvcall_conf_t *conf = (snvcall_conf_t *)confp;
     /* 4 bases ignoring N, -1 reference/consensus base makes 3 */
     long double pvalues[NUM_NONCONS_BASES]; /* pvalues reported back from snpcaller */
     int alt_counts[NUM_NONCONS_BASES]; /* counts for alt bases handed down to snpcaller */
     int alt_raw_counts[NUM_NONCONS_BASES]; /* raw, unfiltered alt-counts */
     int alt_bases[NUM_NONCONS_BASES];/* actual alt bases */
     int got_alt_bases = 0;
     int cons_as_ref = conf->flag & SNVCALL_CONS_AS_REF;

     /* don't call if no coverage or if we don't know what to call
      * against */
     if (p->coverage == 0 || p->cons_base[0] == 'N') {          
          return;
     }
     if (p->coverage < conf->min_cov) {          
          return;
     }
     if (! conf->dont_skip_n && p->ref_base == 'N') {
          return;
     }
     

     /* call indels
      */
     if (! conf->no_indels) {
          /* FIXME make function */

          /* Report consensus indel
           * FIXME: call other indels/substitutions with respect to consensus indel */
          if (p->cons_base[0] == '+') {
               report_cons_ins(p, conf);
               return;
          } 
          if (p->cons_base[0] == '-') {
               report_cons_del(p, conf);
               return;
          }
          /* Call indels */
          if (p->num_dels || p->num_ins) {         
               /*
               LOG_FIXME("%s:%d (p->num_dels=%d p->del_quals=%d"
                         " p->num_ins=%d p->ins_quals.n=%d\n", 
                         p->target, p->pos+1, p->num_dels, p->del_quals.n,
                         p->num_ins, p->ins_quals.n);
               */

               if (p->num_ins && p->ins_quals.n) {
                    ins_event *it, *it_tmp;
                    HASH_ITER(hh_ins, p->ins_event_counts, it, it_tmp) {
                         plp_to_ins_errprobs(&bi_err_probs, &bi_num_err_probs, 
                              p, conf, it->key);  
                         qsort(bi_err_probs, bi_num_err_probs, sizeof(double), dbl_cmp);
                         if (conf->bonf_dynamic) {
                              conf->bonf_indel += 1;
                         }
                         num_indel_tests += 1;
                         if (call_alt_ins(p, bi_err_probs, bi_num_err_probs, conf, it)) {
                              free(bi_err_probs);
                              return;
                         }
                         free(bi_err_probs);
                    }
               }
               if (p->num_dels && p->del_quals.n) {
                    del_event *it, *it_tmp;
                    HASH_ITER(hh_del, p->del_event_counts, it, it_tmp) {
                         plp_to_del_errprobs(&bd_err_probs, &bd_num_err_probs, 
                              p, conf, it->key);
                         qsort(bd_err_probs, bd_num_err_probs, sizeof(double), dbl_cmp);
                         if (conf->bonf_dynamic) {
                              conf->bonf_indel += 1;
                         }
                         num_indel_tests += 1;
                         if (call_alt_del(p, bd_err_probs, bd_num_err_probs, conf, it)) {
                              free(bd_err_probs);
                              return;
                         }
                         free(bd_err_probs);
                    }
               }
          }
     }

     /* call snvs
      */
     if (! conf->only_indels) {
          /* FIXME make function */

          /* 
           * as snv-ref we always report what the user wanted (given
           * ref or determined cons if corresponding flag was given).
           * we always use the determined cons for calculating a quality though.
           *
           *  Example with given ref A and cons/alt C
           *  Alt bases are all non-cons bases: A(!), T, G
           * 
           *  if cons_as_ref:
           *   C>A tested as A vs C
           *   C>T tested as T vs C
           *   C>G tested as G vs C
           *   (no CONSVARs possible by definition)
           *    
           *  else (ref):
           *   A>C tested as NIL (CONSVAR A vs C)
           *   A>T tested as T vs C
           *   A>G tested as G vs C
           *   (A>A obviously not possible)
           *
           * First branch is "default" i.e. doesn't need exception checking
           */

          /* CONSVAR, i.e. the consensus determined here is different from
           * the reference coming from a fasta file
           */
          if (p->ref_base != p->cons_base[0] && !cons_as_ref) {/* && p->ref_base != 'N') {*/
               report_cons_sub(p, conf);
          }

          /* no point continuing we're calling against a reference (not
           * consensus) and ref nt is an n. only consvars should be
           * predicted then (see above) */
          if (p->ref_base == 'N' && !cons_as_ref) {
               return;
          }

          plp_to_errprobs(&bc_err_probs, &bc_num_err_probs, 
                          alt_bases, alt_counts, alt_raw_counts,
                          p, conf);

     #if 0
          for (i=0; i<NUM_NONCONS_BASES; i++) {
               LOG_FIXME("NUM_NONCONS_BASES=%d alt_counts=%d alt_raw_counts=%d\n", i, alt_counts[i], alt_raw_counts[i]);
          }
     #endif

          for (i=0; i<NUM_NONCONS_BASES; i++) {
               if (alt_counts[i]) {
                    got_alt_bases = 1;
                    break;
               }
          }
          if (! got_alt_bases) {
               LOG_DEBUG("%s %d: only cons bases left after filtering.\n", 
                         p->target, p->pos+1);
               /* ...and CONSVAR already reported */
               free(bc_err_probs);
               return;
          }

          /* sorting in ascending order should in theory be numerically
           * more stable and also make snpcaller faster */
          qsort(bc_err_probs, bc_num_err_probs, sizeof(double), dbl_cmp);

     #ifdef TRACE
          {
               int i=0;
               for (i=0; i<bc_num_err_probs; i++) {
                    LOG_FATAL("after sorting i=%d err_prob=%g\n", i, bc_err_probs[i]);
               }
          }
     #endif
          if (conf->bonf_dynamic) {
               if (1 == conf->bonf_sub) {
                    conf->bonf_sub = NUM_NONCONS_BASES; /* otherwise we start with 1+NUM_NONCONS_BASES */
               } else {
                    conf->bonf_sub += NUM_NONCONS_BASES; /* will do one test per non-cons nuc */
               }
          }
          num_snv_tests += NUM_NONCONS_BASES;

          LOG_DEBUG("%s %d: passing down %d quals with noncons_counts"
                    " (%d, %d, %d) to snpcaller(num_snv_tests=%lld conf->bonf=%lld, conf->sig=%f)\n", p->target, p->pos+1,
                    bc_num_err_probs, alt_counts[0], alt_counts[1], alt_counts[2], num_snv_tests, conf->bonf_sub, conf->sig);

          if (snpcaller(pvalues, bc_err_probs, bc_num_err_probs, 
                       alt_counts, conf->bonf_sub, conf->sig)) {
               fprintf(stderr, "FATAL: snpcaller() failed at %s:%s():%d\n",
                       __FILE__, __FUNCTION__, __LINE__);
               free(bc_err_probs);
               return;
          }

          /* for all alt-bases, i.e. non-cons bases (which might include
           * the ref-base!) */
          for (i=0; i<NUM_NONCONS_BASES; i++) {
               int alt_base = alt_bases[i];
               int alt_count = alt_counts[i];
               int alt_raw_count = alt_raw_counts[i];
               long double pvalue = pvalues[i];
               int reported_snv_ref = cons_as_ref ? p->cons_base[0] : p->ref_base;

               if (alt_base==reported_snv_ref) { /* p->ref_base && !cons_as_ref) {*/
                    /* self comparison */
     #if DEBUG
                    LOG_DEBUG("%s\n", "continue because self comparison")
     #endif
                    continue;
               }
               if (p->ref_base==alt_base && !cons_as_ref) {
     #if DEBUG
                    LOG_DEBUG("%s\n", "continue because: p->ref_base==alt_base && !cons_as_ref")
     #endif
                    continue;
               }

               if (pvalue * (double)conf->bonf_sub < conf->sig) {
                    const int is_indel = 0;
                    const int is_consvar = 0;
                    float af = alt_raw_count/(float)p->coverage;

                    char report_ref[2];
                    char report_alt[2];
                    report_ref[0] = reported_snv_ref;
                    report_alt[0] = alt_base;
                    report_ref[1] = report_alt[1] = '\0';

                    int ref_nt4;
                    int alt_nt4;
                    ref_nt4 = bam_nt4_table[(int)report_ref[0]];
                    alt_nt4 = bam_nt4_table[(int)report_alt[0]];

                    dp4_counts_t dp4;
                    dp4.ref_fw = p->fw_counts[ref_nt4];
                    dp4.ref_rv = p->rv_counts[ref_nt4];
                    dp4.alt_fw = p->fw_counts[alt_nt4];
                    dp4.alt_rv = p->rv_counts[alt_nt4];

                    report_var(& conf->vcf_out, p, report_ref, report_alt, 
                               af, PROB_TO_PHREDQUAL(pvalue), 
                               is_indel, is_consvar, &dp4);
                    LOG_DEBUG("low freq snp: %s %d %c>%c pv-prob:%Lg;pv-qual:%d"
                              " counts-raw:%d/%d=%.6f counts-filt:%d/%d=%.6f\n",
                              p->target, p->pos+1, p->cons_base[0], alt_base,
                              pvalue, PROB_TO_PHREDQUAL(pvalue),
                              /* counts-raw */ alt_raw_count, p->coverage, alt_raw_count/(float)p->coverage,
                              /* counts-filt */ alt_count, bc_num_err_probs, alt_count/(float)bc_num_err_probs);
               }
     #if 0
               else {
                    LOG_DEBUG("non sig: pvalue=%Lg * (double)conf->bonf=%lld < conf->sig=%f\n", pvalue, conf->bonf, conf->sig);
               }
     #endif
          }
          free(bc_err_probs);
     }

}
/* call_vars() */



static void
usage(const mplp_conf_t *mplp_conf, const snvcall_conf_t *snvcall_conf)
{
     fprintf(stderr, "Usage: %s [options] in.bam\n\n", MYNAME);
     fprintf(stderr, "Options:\n");
     fprintf(stderr, "- Regions\n");                                        
     fprintf(stderr, "       -r | --region STR            Region in which pileup should be generated [null]\n");
     fprintf(stderr, "       -l | --bed FILE              List of positions (chr pos) or regions (BED) [null]\n");

     fprintf(stderr, "- Reference\n");                                               
     fprintf(stderr, "       -f | --ref FILE            Indexed reference fasta file (gzip supported) [null]\n");
     fprintf(stderr, "            --cons-as-ref           Use consensus base as ref, i.e. ignore base given in reffa (reffa still used for BAQ, if enabled)\n");

     fprintf(stderr, "- Output\n");                                                
     fprintf(stderr, "       -o | --out FILE              Vcf output file [- = stdout]\n");

     fprintf(stderr, "- Base-call quality\n");                      
     fprintf(stderr, "       -q | --min-bq INT            Skip any base with baseQ smaller than INT [%d]\n", snvcall_conf->min_bq);
     fprintf(stderr, "       -Q | --min-alt-bq INT        Skip alternate bases with baseQ smaller than INT [%d]\n", snvcall_conf->min_alt_bq);
     fprintf(stderr, "       -R | --def-alt-bq INT        Overwrite baseQs of alternate bases (that passed bq filter) with this value (-1: use median ref-bq; 0: keep) [%d]\n", snvcall_conf->def_alt_bq);

     fprintf(stderr, "       -j | --min-jq INT            Skip any base with joinedQ smaller than INT [%d]\n", snvcall_conf->min_jq);
     fprintf(stderr, "       -J | --min-alt-jq INT        Skip alternate bases with joinedQ smaller than INT [%d]\n", snvcall_conf->min_alt_jq);
     fprintf(stderr, "       -K | --def-alt-jq INT        Overwrite joinedQs of alternate bases (that passed jq filter) with this value (-1: use median ref-bq; 0: keep) [%d]\n", snvcall_conf->def_alt_jq);

     fprintf(stderr, "- Base-alignment quality (BAQ)\n");                      
     fprintf(stderr, "       -B | --no-baq                Disable use of base-alignment quality (BAQ)\n");
     fprintf(stderr, "       -D | --del-baq               Delete pre-existing BAQ values, i.e. compute even if already present in BAM\n");
     fprintf(stderr, "       -e | --no-ext-baq            Use 'normal' BAQ (samtools default) instead of extended BAQ (both computed on the fly if not already present in %s tag)\n", BAQ_TAG);

     fprintf(stderr, "- Mapping quality\n");                                
     fprintf(stderr, "       -m | --min-mq INT            Skip reads with mapping quality smaller than INT [%d]\n", mplp_conf->min_mq);
     fprintf(stderr, "       -M | --max-mq INT            Cap mapping quality at INT [%d]\n", mplp_conf->max_mq);
     fprintf(stderr, "       -N | --no-mq                 Don't merge mapping quality in LoFreq's model\n");

     fprintf(stderr, "- Indels\n");                                
     fprintf(stderr, "            --only-indels           Only call indels; no SNVs\n");
     fprintf(stderr, "            --no-indels             Don't call indels; only SNVs\n");

#ifdef USE_SOURCEQUAL                                     
     fprintf(stderr, "- Source quality\n");                                
     fprintf(stderr, "       -s | --src-qual           Enable computation of source quality for reads\n");
     fprintf(stderr, "       -S | --ign-vcf FILE          Ignore variants in this vcf file for source quality computation\n"),
     fprintf(stderr, "       -T | --def-nm-q INT          If >= 0, then replace non-match base qualities with this default value [%d]\n", mplp_conf->def_nm_q);
#endif                                                    
     fprintf(stderr, "- P-Values\n");                                          
     fprintf(stderr, "       -a | --sig                   P-Value cutoff / significance level [%f]\n", snvcall_conf->sig);
     fprintf(stderr, "       -b | --bonf                  Bonferroni factor. 'dynamic' (increase per actually performed test) or INT ['dynamic']\n");
     fprintf(stderr, "- Misc.\n");                                           
#ifdef USE_ALNERRPROF
     fprintf(stderr, "       -A | --map-prof FILE         Mapping error profile (produced with bamstats)\n");
#endif
     fprintf(stderr, "       -C | --min-cov INT           Test only positions having at least this coverage [%d]\n", snvcall_conf->min_cov);
     fprintf(stderr, "            --illumina-1.3          Assume the quality is Illumina-1.3-1.7/ASCII+64 encoded\n");
     fprintf(stderr, "            --dont-skip-n           Don't skip positions where refbase is N (will try to predict CONSVARs (only) at those positions)\n");
     fprintf(stderr, "            --use-orphan            Count anomalous read pairs (i.e. where mate is not aligned properly)\n");
     fprintf(stderr, "            --plp-summary-only      No snv-calling: just output pileup summary per column\n");
     fprintf(stderr, "            --no-default-filter     Don't run default 'lofreq filter' automatically after calling variants\n");
     fprintf(stderr, "            --verbose               Be verbose\n");
     fprintf(stderr, "            --debug                 Enable debugging\n");
}
/* usage() */


int 
main_call(int argc, char *argv[])
{
     /* based on bam_mpileup() */
     int c, i;
     static int use_orphan = 0;
     static int cons_as_ref = 0;
     static int only_indels = 0;
     static int no_indels = 0;

     static int plp_summary_only = 0;
     static int no_default_filter = 0;
     static int dont_skip_n = 0;
     static int illumina_1_3 = 0;
     char *bam_file = NULL;
     char *bed_file = NULL;
     char *vcf_out = NULL; /* == - == stdout */
     char vcf_tmp_template[] = "/tmp/lofreq2-call-dyn-bonf.XXXXXX";
     char *vcf_tmp_out = NULL; /* write to this file first, then filter */
     mplp_conf_t mplp_conf;
     snvcall_conf_t snvcall_conf;
     /*void (*plp_proc_func)(const plp_col_t*, const snvcall_conf_t*);*/
     void (*plp_proc_func)(const plp_col_t*, void*);
     int rc = 0;
     char *ign_vcf = NULL;



/* FIXME add sens test:
construct p such with
quality_range = [20, 25, 30, 35, 40]
coverage_range = [10, 50, 100, 500, 1000, 5000, 10000]
refbase = 'A'
snpbase = 'C'
for cov in coverage_range:
    for q in quality_range:
        num_noncons = 1
        while True:
            void call_snvs(const plp_col_t *p, &snvcall_conf);
            count snvs in output
            if len(snps):
                print num_noncons
                break
            num_noncons += 1
            if num_noncons == cov:
                break
*/


     for (i=0; i<argc; i++) {
          LOG_DEBUG("arg %d: %s\n", i, argv[i]);
     }

     /* default pileup options */
     init_mplp_conf(& mplp_conf);

     /* default snvcall options */
     init_snvcall_conf(& snvcall_conf);
    
    /* keep in sync with long_opts_str and usage 
     *
     * getopt is a pain in the whole when it comes to syncing of long
     * and short args and usage. check out gopt, libcfu...
     */
    while (1) {
         static struct option long_opts[] = {
              /* see usage sync */
              {"region", required_argument, NULL, 'r'},
              {"bed", required_argument, NULL, 'l'}, /* changes here must be reflected in pseudo_parallel code as well */
              
              {"ref", required_argument, NULL, 'f'},
              {"cons-as-ref", no_argument, &cons_as_ref, 1},
              {"no-indels", no_argument, &no_indels, 1},
              {"only-indels", no_argument, &only_indels, 1},

              {"out", required_argument, NULL, 'o'}, /* NOTE changes here must be reflected in pseudo_parallel code as well */

              {"min-bq", required_argument, NULL, 'q'},
              {"min-alt-bq", required_argument, NULL, 'Q'},
              {"def-alt-bq", required_argument, NULL, 'R'},

              {"min-jq", required_argument, NULL, 'j'},
              {"min-alt-jq", required_argument, NULL, 'J'},
              {"def-alt-jq", required_argument, NULL, 'K'},

              {"no-baq", no_argument, NULL, 'B'},
              {"del-baq", no_argument, NULL, 'D'},
              {"no-ext-baq", no_argument, NULL, 'e'},
                  
              {"min-mq", required_argument, NULL, 'm'},
              {"max-mq", required_argument, NULL, 'M'},
              {"no-mq", no_argument, NULL, 'N'},
#ifdef USE_SOURCEQUAL
              {"src-qual", no_argument, NULL, 's'},
              {"ign-vcf", required_argument, NULL, 'S'},
              {"def-nm-q", required_argument, NULL, 'R'},
#endif
              {"sig", required_argument, NULL, 'a'},
              {"bonf", required_argument, NULL, 'b'}, /* NOTE changes here must be reflected in pseudo_parallel code as well */

#ifdef USE_ALNERRPROF
              {"map-prof", required_argument, NULL, 'A'},
#endif                   
              {"min-cov", required_argument, NULL, 'C'},
              /*{"maxdepth", required_argument, NULL, 'd'},*/

              {"dont-skip-n", required_argument, &dont_skip_n, 1},
              {"illumina-1.3", no_argument, &illumina_1_3, 1},
              {"use-orphan", no_argument, &use_orphan, 1},
              {"plp-summary-only", no_argument, &plp_summary_only, 1},
              {"no-default-filter", no_argument, &no_default_filter, 1},
              {"verbose", no_argument, &verbose, 1},
              {"debug", no_argument, &debug, 1},
              {"help", no_argument, NULL, 'h'},

              {0, 0, 0, 0} /* sentinel */
         };

         /* keep in sync with long_opts and usage */
         static const char *long_opts_str = "r:l:f:o:q:Q:R:j:J:K:BDem:M:NsS:T:a:b:A:C:h"; 
         
         /* getopt_long stores the option index here. */
         int long_opts_index = 0;
         c = getopt_long(argc-1, argv+1, /* skipping 'lofreq', just leaving 'command', i.e. call */
                         long_opts_str, long_opts, & long_opts_index);
         if (c == -1) {
              break;
         }

         switch (c) {
         /* see usage sync */
         case 'r': 
              mplp_conf.reg = strdup(optarg); 
              /* FIXME you can enter lots of invalid stuff and libbam
               * won't complain. add checks here or late */
              break;

         case 'l': 
              bed_file = strdup(optarg);
              break;

         case 'f':
              if (! file_exists(optarg)) {
                   LOG_FATAL("Reference fasta file '%s' does not exist. Exiting...\n", optarg);
                   return 1;
              }
              mplp_conf.fa = strdup(optarg);
              mplp_conf.fai = fai_load(optarg);
              if (mplp_conf.fai == 0)  {
                   free(mplp_conf.fa);
                   return 1;
              }
              warn_old_fai(mplp_conf.fa);
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

         case 'q': 
              snvcall_conf.min_bq = atoi(optarg); 
              break;

         case 'Q': 
              snvcall_conf.min_alt_bq = atoi(optarg); 
              break;

         case 'R':
              snvcall_conf.def_alt_bq = atoi(optarg); 
              break;

         case 'j': 
              snvcall_conf.min_jq = atoi(optarg); 
              break;

         case 'J': 
              snvcall_conf.min_alt_jq = atoi(optarg); 
              break;

         case 'K':
              snvcall_conf.def_alt_jq = atoi(optarg);
              if (-1 == snvcall_conf.def_alt_jq) {
                   LOG_FATAL("%s\n", "Sorry, use of median ref JQ implemented yet");/* FIXME */
                   exit(1);
              }
              break;

         case 'B': 
              mplp_conf.flag &= ~MPLP_BAQ; 
              snvcall_conf.flag &= ~SNVCALL_USE_BAQ; 
              break;

         case 'D': 
              mplp_conf.flag |= MPLP_REDO_BAQ;
              break;

         case 'e': 
              mplp_conf.flag &= ~MPLP_EXT_BAQ;
              break;

         case 'm': 
              mplp_conf.min_mq = atoi(optarg); 
              break;

         case 'M': 
              mplp_conf.max_mq = atoi(optarg); 
              break;

         case 'N': 
              snvcall_conf.flag &= ~SNVCALL_USE_MQ; 
              break;

#ifdef USE_SOURCEQUAL
         case 's': 
              mplp_conf.flag |= MPLP_USE_SQ;
              snvcall_conf.flag |= SNVCALL_USE_SQ; 
              break;

         case 'S': 
              ign_vcf = strdup(optarg);
              break;

         case 'T': 
              mplp_conf.def_nm_q = atoi(optarg);
              break;
#else
         case 's': 
              LOG_FATAL("%s\n", "source-qual was disabled at compile time"); 
              return 1;
         case 'S': 
              LOG_FATAL("%s\n", "source-qual was disabled at compile time"); 
              return 1;
         case 'R': 
              LOG_FATAL("%s\n", "source-qual was disabled at compile time"); 
              return 1;
#endif

         case 'a': 
              snvcall_conf.sig = strtof(optarg, (char **)NULL); /* atof */
              if (0==snvcall_conf.sig) {
                   LOG_FATAL("%s\n", "Couldn't parse sign-threshold"); 
                   return 1;
              }
              break;
         case 'b': 
              if (0 == strncmp(optarg, "dynamic", 7)) {
                   snvcall_conf.bonf_dynamic = 1;
                   
              } else {
                   snvcall_conf.bonf_dynamic = 0;

                   snvcall_conf.bonf_sub = strtoll(optarg, (char **)NULL, 10); /* atol */ 
                   if (1>snvcall_conf.bonf_sub) {
                        LOG_FATAL("%s\n", "Couldn't parse Bonferroni factor"); 
                        return 1;
                   }
              }
              break;


#ifdef USE_ALNERRPROF
         case 'A':
              mplp_conf.alnerrprof_file = strdup(optarg);
              break;
#endif

         case 'C': 
              snvcall_conf.min_cov = atoi(optarg); 
              break;

         case 'h': 
              usage(& mplp_conf, & snvcall_conf); 
              return 0; /* WARN: not printing defaults if some args where parsed */

         case '?': 
              LOG_FATAL("%s\n", "unrecognized arguments found. Exiting...\n"); 
              free(bed_file);
              free(vcf_out);
              return 1;
#if 0
         case 0:
              fprintf(stderr, "ERROR: long opt (%s) not mapping to short option."
                      " Exiting...\n", long_opts[long_opts_index].name); 
              return 1;
#endif
         default:
              break;
         }
    }


    snvcall_conf.dont_skip_n = dont_skip_n;
    snvcall_conf.no_indels = no_indels;
    snvcall_conf.only_indels = only_indels;

    if (illumina_1_3) {
         mplp_conf.flag |= MPLP_ILLUMINA13; 
    }

    if (use_orphan) {
         mplp_conf.flag &= ~MPLP_NO_ORPHAN;
    }

    if (cons_as_ref) {
         snvcall_conf.flag |= SNVCALL_CONS_AS_REF;
    }

    if (no_indels && only_indels) {
         LOG_FATAL("%s\n", "Invalid user request to predict no-indels *and* only-indels!? Exiting...\n"); 
         return -1;
    }
    
    if (argc == 2) {
        fprintf(stderr, "\n");
        usage(& mplp_conf, & snvcall_conf);
        return 1;
    }

   /* get bam file argument
    */
    if (1 != argc - optind - 1) {
        fprintf(stderr, "Need exactly one BAM file as last argument\n");
        return 1;
    }
    bam_file = (argv + optind + 1)[0];
        if (0 == strcmp(bam_file, "-")) {
         if (mplp_conf.reg) {
              LOG_FATAL("%s\n", "Need index if region was given and"
                        " index file can't be provided when using stdin mode.");
              return 1;
         }
    } else {
         if (! file_exists(bam_file)) {
              LOG_FATAL("BAM file %s does not exist. Exiting...\n", bam_file);
              return 1;
         }
    }


    /* FIXME: implement function for checking user arg logic */
    if (mplp_conf.min_mq > mplp_conf.max_mq) {
         LOG_FATAL("Minimum mapping quality (%d) larger than maximum mapping quality (%d)\n",
                   mplp_conf.min_mq, mplp_conf.max_mq); 
         return 1;
    }
    if (snvcall_conf.min_bq > snvcall_conf.min_alt_bq) {
         LOG_FATAL("Minimum base-call quality for all bases (%d) larger than minimum base-call quality for alternate bases (%d)\n",
                   snvcall_conf.min_bq, snvcall_conf.min_alt_bq); 
         return 1;
    }
    if (mplp_conf.flag & MPLP_BAQ && ! mplp_conf.fa && ! plp_summary_only) {
         LOG_FATAL("%s\n", "Can't compute BAQ with no reference...\n"); 
         return 1;
    }
    if ( ! (snvcall_conf.flag & SNVCALL_CONS_AS_REF) && ! mplp_conf.fa && ! plp_summary_only) {
         LOG_FATAL("%s\n", "Need a reference when not calling in consensus mode...\n"); 
         return 1;
    }


    if (! plp_summary_only & ! mplp_conf.fa) {
         LOG_WARN("%s\n", "Calling SNVs without reference\n"); 
    }

  
    /* if we don't apply a default filter and bonf is not dynamic then
     * we can directly write to requested output file. otherwise we
     * use a tmp file that gets filtered.
     */
    if (no_default_filter && ! snvcall_conf.bonf_dynamic) {
         if (NULL == vcf_out || 0 == strcmp(vcf_out, "-")) {
              if (vcf_file_open(& snvcall_conf.vcf_out, "-", 
                                0, 'w')) {
                   LOG_ERROR("%s\n", "Couldn't open stdout");
                   return 1;
              }
         } else {
              if (vcf_file_open(& snvcall_conf.vcf_out, vcf_out,
                                HAS_GZIP_EXT(vcf_out), 'w')) {
                   LOG_ERROR("Couldn't open %s\n", vcf_out);
                   return 1;
              }
         }
    } else {
         vcf_tmp_out = strdup(mktemp(vcf_tmp_template));
         if (NULL == vcf_tmp_out) {
              LOG_FATAL("%s\n", "Couldn't create temporary vcf file");
              return 1;
         }
         if (vcf_file_open(& snvcall_conf.vcf_out, vcf_tmp_out,
                           HAS_GZIP_EXT(vcf_tmp_out), 'w')) {
              LOG_ERROR("Couldn't open %s\n", vcf_tmp_out);
              return 1;
         }
    }


    /* save command-line for later reference */
    mplp_conf.cmdline[0] = '\0';
    for (i=0; i<argc; i++) {
         strncat(mplp_conf.cmdline, argv[i], 
                 sizeof(mplp_conf.cmdline)-strlen(mplp_conf.cmdline)-2);
         strcat(mplp_conf.cmdline, " ");
    }

    if (bed_file) {
         mplp_conf.bed = bed_read(bed_file);
         if (! mplp_conf.bed) {
              LOG_ERROR("Couldn't read %s\n", bed_file);
              return 1;
         }
    }

#ifdef USE_SOURCEQUAL
    if (ign_vcf) {
         /* FIXME ignore variants outside given region (on top of bed as well) */
         if (source_qual_load_ign_vcf(ign_vcf, mplp_conf.bed)) {
              LOG_FATAL("Loading of ignore positions from %s failed.", optarg);
              return 1;
         }
         free(ign_vcf);
    }
#endif

    if (debug) {
         dump_mplp_conf(& mplp_conf, stderr);
         dump_snvcall_conf(& snvcall_conf, stderr);
    }

#if 0
    LOG_FIXME("%s\n", "Loading hardcoded vcf file to ignore for source_qual()");
    if (source_qual_load_ign_vcf("./tests/schmock.vcf")) {
         LOG_FATAL("Loading of ignore positions from %s failed.", "FIXME");
         return 1;
    }
#endif

    if (plp_summary_only) {
         plp_proc_func = &plp_summary;

    } else {
         /* or use PACKAGE_STRING */
         vcf_write_new_header(& snvcall_conf.vcf_out,
                              mplp_conf.cmdline, mplp_conf.fa);
         plp_proc_func = &call_vars;
    }

    rc = mpileup(&mplp_conf, plp_proc_func, (void*)&snvcall_conf,
                 1, (const char **) argv + optind + 1);
    if (rc) {
         return rc;
    }

    vcf_file_close(& snvcall_conf.vcf_out);

    /* snv calling completed. now filter according to the following rules:
     *  1. no_default_filter and ! dyn
     *     just print
     *  2 filter with
     *     - no_default_filter, if set
     *     - filter snvphred according to bonf, if dynamic
     */
    if (plp_summary_only) {
         LOG_VERBOSE("%s\n", "No filtering needed: didn't run in SNV calling mode");

    } else if (no_default_filter && ! snvcall_conf.bonf_dynamic) {
         /* vcf file needs no filtering and was already printed to
          * final destination. already taken care of above. */
         LOG_VERBOSE("%s\n", "No filtering needed or requested: variants already written to final destination");

    } else {
         char cmd[BUF_SIZE];
         int len; 

         snprintf(cmd, BUF_SIZE, 
                  "lofreq filter --only-passed -i %s -o %s",
                  vcf_tmp_out, NULL==vcf_out ? "-" : vcf_out);
         len = strlen(cmd);

         if (no_default_filter) {
              len += sprintf(cmd+len, " %s", "--no-defaults");
         }

         if (snvcall_conf.bonf_dynamic) {
              len += sprintf(cmd+len,/* appending to str with format. see http://stackoverflow.com/questions/14023024/strcat-for-formatted-strings */
                             " --snvqual-thresh %d --indelqual-thresh %d", 
                             snvcall_conf.bonf_sub ? PROB_TO_PHREDQUAL(snvcall_conf.sig/snvcall_conf.bonf_sub) : INT_MAX,
                             snvcall_conf.bonf_indel ? PROB_TO_PHREDQUAL(snvcall_conf.sig/snvcall_conf.bonf_indel) : INT_MAX);
         } else {
              LOG_VERBOSE("%s\n", "No SNV/indel-quality filtering needed (already applied during call since bonf was fixed)");
         }

         LOG_VERBOSE("Executing %s\n", cmd);
         if (0 != (rc = system(cmd))) {
              LOG_ERROR("The following command failed: %s\n", cmd);
              rc = 1;
              
         } else {
              /*if (! debug)*/
              (void) unlink(vcf_tmp_out);
         }
    }

    if (! plp_summary_only) {
         /* output some stats. number of tests performed need for
          * multiple testing correction. line will be parse by
          * downstream script e.g. lofreq_somatic, so be careful when
          * changing the format */
         int org_verbose = verbose;
         verbose = 1;
         /* lofreq2_call_parallel.py and used by lofreq2_somatic.py */
         LOG_VERBOSE("Number of substitution tests performed: %lld\n", num_snv_tests);
         LOG_VERBOSE("Number of indel tests performed: %lld\n", num_indel_tests);
         verbose = org_verbose;
    }

#ifdef USE_SOURCEQUAL
    source_qual_free_ign_vars();
#endif

    free(vcf_tmp_out);
    free(vcf_out);
    free(mplp_conf.alnerrprof_file);
    free(mplp_conf.reg); 
    free(mplp_conf.fa);
    if (mplp_conf.fai) {
         fai_destroy(mplp_conf.fai);
    }
    free(bed_file);
    if (mplp_conf.bed) {
         bed_destroy(mplp_conf.bed);
    }
    
    if (0==rc) {
         LOG_VERBOSE("%s\n", "Successful exit.");
    }

    return rc;
}
/* main_call */
