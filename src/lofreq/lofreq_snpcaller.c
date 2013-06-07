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


#define SNVCALL_USE_MQ      0x10
#define SNVCALL_USE_SQ      0x20
#define SNVCALL_CONS_AS_REF 0x40

#define BUF_SIZE 1<<16


/* scale 0-60 to from 0-254 and shrink
 * Y = 254.0/60.0 * MQ * (MQ**X)/(60**X)
 *
 * if 20 should stay 20
 * 20 = 254/60.0 * 20 * (20**X)/(60**X) 
 * 60/254.0 = (20**X)/(60**X)
 * (20/60.0)**X = 60/254.0
 * since a**x = y equals log_a(y) = x
 * x = log_a(60/254.0); a=20/60.0;
 * x = 1.3134658329243962 
 */
#undef SCALE_MQ
#define SCALE_MQ_FAC  1.3134658329243962

/* filled in missing values with the min of the two neighbouring values */
#define TRUE_MQ_BWA_HG19_EXOME_2X100_SIMUL
#undef TRUE_MQ_BWA_HG19_EXOME_2X100_SIMUL
#ifdef TRUE_MQ_BWA_HG19_EXOME_2X100_SIMUL
const int MQ_TRANS_TABLE[61] = {
1,
1,
3,
4,
5,
5,
8,
9,
4,
8,
14,
17,
22,
25,
25,
29,
32,
33,
34,
34, /* NA */
34,
34, /* NA */
34, /* NA */
34,
34, /* NA */
34, /* NA */
34, /* NA */
34, /* NA */
34, /* NA */
41,
41, /* NA */
41, /* NA */
41, /* NA */
41, /* NA */
41, /* NA */
41, /* NA */
41, /* NA */
41, /* NA */
50,
46, /* NA */
46, /* NA */
46, /* NA */
46, /* NA */
46, /* NA */
46, /* NA */
46, /* NA */
46, /* NA */
46, /* NA */
46, /* NA */
46, /* NA */
46, /* NA */
46,
46, /* NA */
46, /* NA */
54,
37,
37, /* NA */
45,
45, /* NA */
45, /* NA */
67};
#endif


typedef struct {
     int min_altbq, def_altbq;
     int bonf_dynamic; /* boolean: incr bonf as we go along. eventual
                        * filtering of all has to be done by
                        * caller! */
     int min_cov;
     long long int bonf; /* warning: changed dynamically ! */
     float sig;
     FILE *out;
     int flag;
} snvcall_conf_t;


void
report_var(FILE *stream, const plp_col_t *p, const char ref, 
           const char alt, const float af, const int qual,
           const int is_indel, const int is_consvar)
{
     var_t *var;
     dp4_counts_t dp4;
     double sb_prob, sb_left_pv, sb_right_pv, sb_two_pv;
     int sb_qual;
     
     vcf_new_var(&var);
     var->chrom = strdup(p->target);
     var->pos = p->pos;
     /* var->id = NA */
     var->ref = ref;
     var->alt = alt;
     if (qual>-1) {
          var->qual = qual;
     }
     /* var->filter = NA */ 
   
     dp4.ref_fw = p->fw_counts[bam_nt4_table[(int)ref]];
     dp4.ref_rv = p->rv_counts[bam_nt4_table[(int)ref]];
     dp4.alt_fw = p->fw_counts[bam_nt4_table[(int)alt]];
     dp4.alt_rv = p->rv_counts[bam_nt4_table[(int)alt]];

     sb_prob = kt_fisher_exact(dp4.ref_fw, dp4.ref_rv, 
                               dp4.alt_fw, dp4.alt_rv,
                               &sb_left_pv, &sb_right_pv, &sb_two_pv);
     sb_qual = PROB_TO_PHREDQUAL(sb_two_pv);

     vcf_var_sprintf_info(var, &p->coverage, &af, &sb_qual,
                          &dp4, is_indel, is_consvar);
     vcf_write_var(stream, var);
     vcf_free_var(&var);
}
/* report_var() */



/* "Merge" MQ and BQ if requested using the following equation:
 *  P_jq = P_mq * + (1-P_mq) P_bq.
 *  If MQ==255 (i.e. not available) return BQ
 */
double
merge_baseq_and_mapq(const int bq, const int mq)
{
     double mp, bp, jp; /* corresponding probs */

     bp = PHREDQUAL_TO_PROB(bq);
     if (255 == mq) {
          return bp;
     }

#ifdef SCALE_MQ
     mp = PHREDQUAL_TO_PROB(254/60.0*mq * pow(mq, SCALE_MQ_FAC)/pow(60, SCALE_MQ_FAC));
#else
     mp = PHREDQUAL_TO_PROB(mq);
#endif
      
#ifdef TRUE_MQ_BWA_HG19_EXOME_2X100_SIMUL
     assert(mq <= 60);
     mp = PHREDQUAL_TO_PROB(MQ_TRANS_TABLE[mq]); 
#endif
     /* No need to do computation in phred-space as
      * numbers won't get small enough.
      */

     /* note: merging Q1 with anything else will result in Q0. */
     jp = mp + (1.0 - mp) * bp;
#ifdef DEBUG
     LOG_DEBUG("P_M + (1-P_M) P_B:   %g + (1.0 - %g) * %g = %g  ==  Q%d + (1.0 - Q%d) * Q%d  =  Q%d\n",
               mp, mp, bp, jp, mq, mq, bq, jq);
#endif
#if 0
     LOG_DEBUG("BQ %d after merging with MQ %d = %d\n", bq, mq, jq);
#endif

     return jp;
}
/* merge_baseq_and_mapq() */



void
plp_summary(const plp_col_t *plp_col, void* confp) 
{
     FILE* stream = stdout;
     int i;

     fprintf(stream, "%s\t%d\t%c\t%c", plp_col->target, plp_col->pos+1,
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
#if 0
     LOG_FIXME("%s\n", "unfinished");
#endif
}




/* low-freq vars always called against cons_base, which might be
 * different from ref_base. if cons_base != ref_base then it's a
 * cons-var.
 * 
 * Assuming conf->min_bq and read-level filtering was already done
 * upstream. altbase mangling happens here however.
 * 
 */
void
call_snvs(const plp_col_t *p, void *confp)
{
     double *err_probs; /* error probs (qualities) passed down to snpcaller */
     int num_err_probs; /* #elements in err_probs */
     int i, j;

     snvcall_conf_t *conf = (snvcall_conf_t *)confp;
     /* 4 bases ignoring N, -1 reference/consensus base makes 3 */
     double pvalues[3]; /* pvalues reported back from snpcaller */
     int alt_counts[3]; /* counts for alt bases handed down to snpcaller */
     int alt_raw_counts[3]; /* raw, unfiltered alt-counts */
     int alt_bases[3];/* actual alt bases */
     int alt_idx;
     int got_alt_bases = 0;

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
     int cons_as_ref = conf->flag & SNVCALL_CONS_AS_REF;

     /* don't call if no coverage or if we don't know what to call
      * against */
     if (p->coverage == 0 || p->cons_base == 'N') {          
          return;
     }
     if (p->coverage < conf->min_cov) {          
          return;
     }

     if (conf->bonf_dynamic) {
          conf->bonf += 3; /* will do one test per non-cons nuc */
     }
#ifdef FIXME
     if (p->num_dels || p->num_ins) {
          LOG_FIXME("%s:%d (p->num_dels=%d p->del_quals=%d"
                    " p->num_ins=%d p->ins_quals.n=%d\n", 
                    p->target, p->pos+1, p->num_dels, p->del_quals.n,
                    p->num_ins, p->ins_quals.n);
          if (p->num_dels && p->del_quals.n) {
               LOG_FIXME("Call deletions at %s:%d\n", p->target, p->pos+1);
          }
          if (p->num_ins && p->ins_quals.n) {
               LOG_FIXME("Call insertions at %s:%d\n", p->target, p->pos+1);
          }
     }
#endif
     /* CONSVAR, i.e. the consensus determined here is different from
      * the reference coming from a fasta file 
      */
     if (p->ref_base != p->cons_base && !cons_as_ref && p->ref_base != 'N') {
          const int is_indel = 0;
          const int is_consvar = 1;
          const int qual = -1;
          float af = base_count(p, p->cons_base) / (float)p->coverage;

          report_var(conf->out, p, p->ref_base, p->cons_base,
                     af, qual, is_indel, is_consvar);
          LOG_DEBUG("cons var snp: %s %d %c>%c\n",
                    p->target, p->pos+1, p->ref_base, p->cons_base);          
     }

     if (NULL == (err_probs = malloc(p->coverage * sizeof(double)))) {
          /* coverage = base-count after read level filtering */
          fprintf(stderr, "FATAL: couldn't allocate memory at %s:%s():%d\n",
                  __FILE__, __FUNCTION__, __LINE__);
          free(err_probs);
          return;
     }
    
     num_err_probs = 0;
     alt_idx = -1;
     for (i=0; i<NUM_NT4; i++) {
          int is_alt_base;
          int nt = bam_nt4_rev_table[i];
          if (nt == 'N') {
               continue;
          }

          is_alt_base = 0;
          if (nt != p->cons_base) {
               is_alt_base = 1;

               alt_idx += 1;
               alt_bases[alt_idx] = nt;
               alt_counts[alt_idx] = 0;
               alt_raw_counts[alt_idx] = 0;
          }

          for (j=0; j<p->base_quals[i].n; j++) {
               int bq, mq;
               double final_err_prob; /* == final quality used for snv calling */
#ifdef USE_SOURCEQUAL
               int sq;
#endif
               assert(p->fw_counts[i] + p->rv_counts[i] == p->base_quals[i].n);
               assert(p->base_quals[i].n == p->map_quals[i].n);
               /* FIXME assert(plp_col.map_quals[i].n == plp_col.source_quals[i].n); */
            
               bq = p->base_quals[i].data[j];
               mq = p->map_quals[i].data[j];
               /* FIXME sq = p->source_quals[i].data[j]; */
               
               if (is_alt_base) {
                    alt_raw_counts[alt_idx] += 1;
                    if (bq < conf->min_altbq) {
                         continue; 
                         /* WARNING base counts now invalid. We use
                          * them for freq reporting anyway, otherwise
                          * heterozygous calls are odd */
                    }
                    bq = conf->def_altbq;
                    alt_counts[alt_idx] += 1;
               }

               if ((conf->flag & SNVCALL_USE_MQ)) {
                    final_err_prob = merge_baseq_and_mapq(bq, mq);

               } else {
                    final_err_prob = PHREDQUAL_TO_PROB(bq);
               }

               err_probs[num_err_probs++] = final_err_prob;
          }
     }

     for (i=0; i<3; i++) {
          if (alt_counts[i]) {
               got_alt_bases = 1;
               break;
          }
     }
     if (! got_alt_bases) {
          LOG_DEBUG("%s %d: only cons bases left after filtering.\n", 
                    p->target, p->pos+1);
          /* ...and CONSVAR already reported */
          free(err_probs);
          return;
     }

     /* sorting in ascending order should in theory be numerically
      * more stable and also make snpcaller faster */
     qsort(err_probs, num_err_probs, sizeof(double), dbl_cmp);
     
#ifdef TRACE
     {
          int i=0;
          for (i=0; i<num_err_probs; i++) {
               LOG_FATAL("after sorting i=%d err_prob=%g\n", i, err_probs[i]);
          }

     }
#endif
     LOG_DEBUG("%s %d: passing down %d quals with noncons_counts"
               " (%d, %d, %d) to snpcaller()\n", p->target, p->pos+1,
               num_err_probs, alt_counts[0], alt_counts[1], alt_counts[2]);

     if (snpcaller(pvalues, err_probs, num_err_probs, 
                  alt_counts, conf->bonf, conf->sig)) {
          fprintf(stderr, "FATAL: snpcaller() failed at %s:%s():%d\n",
                  __FILE__, __FUNCTION__, __LINE__);
          free(err_probs);
          return;
     }

     /* for all alt-bases, i.e. non-cons bases (which might include
      * the ref-base!) */
     for (i=0; i<3; i++) {
          int alt_base = alt_bases[i];
          int alt_count = alt_counts[i];
          int alt_raw_count = alt_raw_counts[i];
          double pvalue = pvalues[i];
          int reported_snv_ref = cons_as_ref ? p->cons_base : p->ref_base;

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

          if (pvalue * (double)conf->bonf < conf->sig) {
               const int is_indel = 0;
               const int is_consvar = 0;
               float af = alt_raw_count/(float)p->coverage;

               report_var(conf->out, p, reported_snv_ref, alt_base, 
                          af, PROB_TO_PHREDQUAL(pvalue), 
                          is_indel, is_consvar);
               LOG_DEBUG("low freq snp: %s %d %c>%c pv-prob:%g;pv-qual:%d"
                         " counts-raw:%d/%d=%.6f counts-filt:%d/%d=%.6f\n",
                         p->target, p->pos+1, p->cons_base, alt_base,
                         pvalue, PROB_TO_PHREDQUAL(pvalue),
                         /* counts-raw */ alt_raw_count, p->coverage, alt_raw_count/(float)p->coverage,
                         /* counts-filt */ alt_count, num_err_probs, alt_count/(float)num_err_probs);
          }
#if DEBUG
          else {
               LOG_DEBUG("non sig: pvalue=%g * (double)conf->bonf=%lld < conf->sig=%f\n", pvalue, conf->bonf, conf->sig);
          }
#endif
     }
     free(err_probs);
}
/* call_snvs() */


char *
cigar_from_bam(const bam1_t *b) {
     /* from char *bam_format1_core(const bam_header_t *header, const bam1_t *b, int of) */
     const bam1_core_t *c = &b->core;
     kstring_t str;
     int i;
     str.l = str.m = 0; str.s = 0;
     for (i = 0; i < c->n_cigar; ++i) {
          kputw(bam1_cigar(b)[i]>>BAM_CIGAR_SHIFT, &str);
          kputc("MIDNSHP=X"[bam1_cigar(b)[i]&BAM_CIGAR_MASK], &str);
     }
     return str.s;
}
/* cigar_from_bam() */



/* Count matches and mismatches for an aligned read and also return
 * the corresponding qualities. returns NULL on error or pointer to
 * qual array for n_match and n_mismatch (sum is size). allocated
 * here. user must free
 */
int *
count_matches(int *n_matches, int *n_mismatches,
              const bam1_t *b, const char *ref)
{
     /* modelled after bam.c:bam_calend(), bam_format1_core() and
      * pysam's aligned_pairs 
      */
     uint32_t *cigar = bam1_cigar(b);
     const bam1_core_t *c = &b->core;
     uint32_t pos = c->pos; /* pos on genome */
     uint32_t qpos = 0; /* pos on read/query */
     uint32_t k, i;
     int *quals = NULL;
     int n_quals = 0;
     /* read length */
     int32_t qlen = (int32_t) bam_cigar2qlen(c, bam1_cigar(b));

     *n_matches = 0;
     *n_mismatches = 0;

     if (NULL==ref) {
          return NULL;
     }

     if (NULL == (quals = malloc(qlen * sizeof(int)))) {
          LOG_FATAL("%s\n", "couldn't allocate memory");
          return NULL;
     }

     if (0) {
          fprintf(stderr, "SOURCEQUAL: core.pos %d - calend %d - cigar %s",
                  b->core.pos, bam_calend(&b->core, bam1_cigar(b)), cigar_from_bam(b));
     }
     
     /* loop over cigar to get aligned bases and matches/mismatches
      * and their quals.
      *
      * read: bam_format1_core(NULL, b, BAM_OFDEC);
      */
     for (k=0; k < c->n_cigar; ++k) { /* n_cigar: number of cigar operations */
          int op = cigar[k] & BAM_CIGAR_MASK; /* the cigar operation */
          uint32_t l = cigar[k] >> BAM_CIGAR_SHIFT;
          
          /* following conditionals could be collapsed to much shorter
           * code, but we keep them as they were in pysam's
           * aligned_pairs to make later handling of indels easier
           */
          if (op == BAM_CMATCH) {
               for (i=pos; i<pos+l; i++) {                             
#if 0
                    printf("qpos,i = %d,%d\n", qpos, i);
#endif
                    char ref_nt = ref[i];
                    char read_nt = bam_nt16_rev_table[bam1_seqi(bam1_seq(b), qpos)];
                    int bq = bam1_qual(b)[qpos];
                    
                    if (ref_nt == read_nt) {
                         *n_matches += 1;
                    } else {
                         *n_mismatches += 1;
                    }
                    quals[n_quals++] = bq;

                    qpos += 1;
               }
               pos += l;
               
          } else if (op == BAM_CINS) {
               for (i=pos; i<pos+l; i++) {
#if 0
                    printf("qpos,i = %d,None\n", qpos);
#endif
                    qpos += 1;
               }
               qpos += l;
               
          } else if (op == BAM_CDEL || op == BAM_CREF_SKIP) {
               for (i=pos; i<pos+l; i++) {
#if 0
                    printf("qpos,i = None,%d\n", i);
#endif
               }
               pos += l;
          }
     } /* for k */
     assert(pos == bam_calend(&b->core, bam1_cigar(b))); /* FIXME correct assert? what if clipped? */

     if (0) {
          fprintf(stderr, " - matches %d - mismatches %d\n", *n_matches, *n_mismatches);                                       
     }
     assert(*n_matches + *n_mismatches == n_quals);

     return quals;
}
/* count_matches() */



#if 0
/* Estimate as to how likely it is that this read, given the mapping,
 * comes from this reference genome. P(r not from g|mapping) = 1 - P(r
 * from g). Use qualities of all bases for and poisson-binomial dist
 * (as for core SNV calling). Assumed independence of errors okay: if
 * they are not independent, then the assumption is conservative. Keep
 * all qualities as they are, i.e. don’t replace mismatches with lower
 * values. Rationale: higher SNV quals, means higher chance SNVs are
 * real, therefore higher prob. read does not come from genome. 
 *
 * FIXME: should always ignore heterozygous or known SNV pos!
 *
 * Returns -1 on error. otherwise phred score of source error prob.
 *
 * FIXME: old definition above and below in source
 *
 */
int
source_qual(const bam1_t *b, const char *ref)
{
     double *probvec;
     int src_qual = 255;
     double src_pvalue;
     int *quals;
     int n_matches = 0;
     int n_mismatches = 0;
     int n_quals = 0;

     quals = count_matches(&n_matches, &n_mismatches, b, ref);
     if (NULL == quals) {
          return -1;
     }
     n_quals = n_matches + n_mismatches;

     /* sorting in theory should be numerically more stable and also
      * make snpcallerfaster */
     qsort(quals, n_quals, sizeof(int), int_cmp);
     probvec = poissbin(&src_pvalue, quals,
                        n_quals, n_mismatches, 1.0, 0.05);


     if (src_pvalue>1.0) {/* DBL_MAX is default return value */
          src_pvalue = 1.0;/*-DBL_EPSILON;*/
     }

     LOG_FIXME("src_pvalue = %g. Actual prob = %g\n", 
               src_pvalue, exp(probvec[n_mismatches]));

     /* src_pvalue: what's the prob of seeing n_mismatches or more by
      * chance, given quals? or: how likely is this read from the
      * genome. 1-src_value = prob read is not from genome
      */
     if (0) {
          LOG_FIXME("Orig src_pv = %f", src_pvalue);
     }
     src_pvalue = 1.0-src_pvalue;
     free(probvec);

     src_qual =  PROB_TO_PHREDQUAL(src_pvalue);

     if (0) {
          int i;
          fprintf(stderr, "| src_pv = %f = Q%d for %d/%d mismatches. All quals: ", 
                  src_pvalue, src_qual, n_mismatches, n_quals);
          for (i=0; i<n_quals; i++) {
               fprintf(stderr, " %d", quals[i]);
          }
          fprintf(stderr, "\n");
     }

#if 0
"
PJ = joined Q
PM = map Q
PG = genome Q
PS = source Q


PJ = PM  +  (1-PM) * PG  +  (1-PM) * (1-PG) * PB
# note: niranjan used PS and meant PB I think
# mapping error
# OR
# no mapping error AND genome error
# OR
# no mapping error AND no genome error AND base-error


PJ = PM + (1-PM) * PB
# mapping error OR no mapping error AND base-error
"
#endif
     free(quals);

     return src_qual;
}
/* source_qual() */
#endif



void
dump_snvcall_conf(const snvcall_conf_t *c, FILE *stream) 
{
     fprintf(stream, "snvcall options\n");
     fprintf(stream, "  min_altbq      = %d\n", c->min_altbq);
     fprintf(stream, "  def_altbq      = %d\n", c->def_altbq);
     fprintf(stream, "  min_cov        = %d\n", c->min_cov);
     fprintf(stream, "  bonf           = %lld  (might get recalculated)\n", c->bonf);
     fprintf(stream, "  bonf_dynamic   = %d\n", c->bonf_dynamic);
     fprintf(stream, "  sig            = %f\n", c->sig);
     fprintf(stream, "  out            = %p\n", (void*)c->out);
     fprintf(stream, "  flag & SNVCALL_USE_MQ      = %d\n", c->flag&SNVCALL_USE_MQ?1:0);
#ifdef USE_SOURCEQUAL
     fprintf(stream, "  flag & SNVCALL_USE_SQ      = %d\n", c->flag&SNVCALL_USE_SQ?1:0);
#endif
     fprintf(stream, "  flag & SNVCALL_CONS_AS_REF = %d\n", c->flag&SNVCALL_CONS_AS_REF?1:0);
}


static void
usage(const mplp_conf_t *mplp_conf, const snvcall_conf_t *snvcall_conf)
{
     fprintf(stderr, "Usage: %s call [options] in.bam\n\n", PACKAGE);
     fprintf(stderr, "Options:\n");
     fprintf(stderr, "- Regions\n");                                        
     fprintf(stderr, "       -r | --region STR            region in which pileup should be generated [null]\n");
     fprintf(stderr, "       -l | --bed FILE              list of positions (chr pos) or regions (BED) [null]\n");
     fprintf(stderr, "- Reference\n");                                               
     fprintf(stderr, "       -f | --reffa FILE            faidx indexed reference sequence file [null]\n");
     fprintf(stderr, "       -c | --cons-as-ref           Use consensus base as ref, i.e. ignore base given in reffa (reffa still used for BAQ, if enabled)\n");
     fprintf(stderr, "- Output\n");                                                
     fprintf(stderr, "       -o | --out FILE              vcf output file [- = stdout]\n");
     fprintf(stderr, "- Base-call quality\n");                      
     fprintf(stderr, "       -q | --min-bq INT            skip any base with baseQ smaller than INT [%d]\n", mplp_conf->min_bq);
     fprintf(stderr, "       -Q | --min-altbq INT         skip nonref-bases with baseQ smaller than INT [%d]. Not active if ref is N\n", snvcall_conf->min_altbq);
     fprintf(stderr, "       -a | --def-altbq INT         nonref base qualities will be replace with this value [%d]\n", snvcall_conf->def_altbq);
     /*fprintf(stderr, "       -B | --no-baq                disable BAQ computation\n");*/
     fprintf(stderr, "       -E | --baq                   enable (extended) per-base alignment quality (BAQ) computation\n");
     fprintf(stderr, "- Mapping quality\n");                                
     fprintf(stderr, "       -m | --min-mq INT            skip alignments with mapQ smaller than INT [%d]\n", mplp_conf->min_mq);
     fprintf(stderr, "       -M | --max-mq INT            cap mapping quality at INT [%d]\n", mplp_conf->max_mq);
     fprintf(stderr, "       -J | --no-mq                 don't merge mapQ into baseQ: P_e = P_mq + (1-P_mq) P_bq\n");
#ifdef USE_SOURCEQUAL                                     
     fprintf(stderr, "       -S | --no-sq                 don't merge sourceQ into baseQ\n");
#endif                                                    
     fprintf(stderr, "- P-Values\n");                                          
     fprintf(stderr, "       -s | --sig                   P-Value cutoff / significance level [%f]\n", snvcall_conf->sig);
     fprintf(stderr, "       -b | --bonf                  Bonferroni factor. 'dynamic' (increase per actually performed test), 'auto' (infer from bed-file) or INT ['dynamic']\n");
     fprintf(stderr, "- Misc.\n");                                           
     fprintf(stderr, "       -C | --min-cov INT           Test only positions having at least this coverage [%d]\n", snvcall_conf->min_cov);
     fprintf(stderr, "       -I | --illumina-1.3          assume the quality is Illumina-1.3-1.7/ASCII+64 encoded\n");
     fprintf(stderr, "            --use-orphan            count anomalous read pairs\n");
     fprintf(stderr, "            --plp-summary-only      no snv-calling. just output pileup summary per column\n");
     fprintf(stderr, "            --no-default-filter     Don't apply default filter command after calling variants\n");
     fprintf(stderr, "            --verbose               be verbose\n");
     fprintf(stderr, "            --debug                 enable debugging\n");
}
/* usage() */




int 
main_call(int argc, char *argv[])
{
     /* based on bam_mpileup() */
     int c, i;
     static int use_orphan = 0;
     static int plp_summary_only = 0;
     static int no_default_filter = 0;
     int bonf_auto = 0;
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

#ifdef SCALE_MQ
     LOG_WARN("%s\n", "MQ scaling switched on!");
#endif
#ifdef TRUE_MQ_BWA_HG19_EXOME_2X100_SIMUL
     LOG_WARN("%s\n", "MQ translation switched on!");
#endif
     for (i=0; i<argc; i++) {
          LOG_DEBUG("arg %d: %s\n", i, argv[i]);
     }


     /* default pileup options */
     memset(&mplp_conf, 0, sizeof(mplp_conf_t));
     mplp_conf.max_mq = 255;
     mplp_conf.min_bq = 3;
     mplp_conf.capQ_thres = 0;
     mplp_conf.max_depth = 1000000;
     mplp_conf.flag = MPLP_NO_ORPHAN; /* | MPLP_REALN | MPLP_REDO_BAQ; */
    
     /* default snvcall options */
     memset(&snvcall_conf, 0, sizeof(snvcall_conf_t));
     snvcall_conf.min_altbq = 20;
     snvcall_conf.def_altbq = snvcall_conf.min_altbq;
     snvcall_conf.min_cov = 1;
     snvcall_conf.bonf_dynamic = 1;
     snvcall_conf.bonf = 1;
     snvcall_conf.sig = 0.05;
     snvcall_conf.out = stdout;
     snvcall_conf.flag = SNVCALL_USE_MQ;/* | MPLP_USE_SQ; FIXME */


    /* keep in sync with long_opts_str and usage 
     *
     * getopt is a pain in the whole when it comes to syncing of long
     * and short args and usage. check out libcfu (also has hash
     * functions etc)
     */
    while (1) {
         static struct option long_opts[] = {
              /* see usage sync */
              {"region", required_argument, NULL, 'r'},
              {"bed", required_argument, NULL, 'l'}, /* changes here must be reflected in pseudo_parallel code as well */
              
              {"reffa", required_argument, NULL, 'f'},
              {"cons-as-ref", no_argument, NULL, 'c'},

              {"out", required_argument, NULL, 'o'}, /* NOTE changes here must be reflected in pseudo_parallel code as well */

              {"min-bq", required_argument, NULL, 'q'},
              {"min-altbq", required_argument, NULL, 'Q'},
              {"def-altbq", required_argument, NULL, 'a'},
              /*{"no-baq", no_argument, NULL, 'B'},*/
              {"baq", required_argument, NULL, 'E'},
                   
              {"min-mq", required_argument, NULL, 'm'},
              {"max-mq", required_argument, NULL, 'M'},
              {"no-mq", no_argument, NULL, 'J'},
#ifdef USE_SOURCEQUAL
              {"no-sq", no_argument, NULL, 'S'},
#endif
              {"bonf", required_argument, NULL, 'b'}, /* NOTE changes here must be reflected in pseudo_parallel code as well */
              {"sig", required_argument, NULL, 's'},
                   
              {"min-cov", required_argument, NULL, 'C'},
              /*{"maxdepth", required_argument, NULL, 'd'},*/
              {"illumina-1.3", no_argument, NULL, 'I'},
              {"use-orphan", no_argument, &use_orphan, 1},
              {"plp-summary-only", no_argument, &plp_summary_only, 1},
              {"no-default-filter", no_argument, &no_default_filter, 1},
              {"verbose", no_argument, &verbose, 1},
              {"debug", no_argument, &debug, 1},
              {"help", no_argument, NULL, 'h'},

              {0, 0, 0, 0} /* sentinel */
         };


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

         /* keep in sync with long_opts and usage */
         static const char *long_opts_str = "r:l:f:co:q:Q:a:Em:M:Js:b:C:Ih"; 
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
              break; /* FIXME you can enter lots of invalid stuff and libbam won't complain. add checks here? */

         case 'l': 
              bed_file = strdup(optarg);
              break;

         case 'f':
              mplp_conf.fa = strdup(optarg);
              mplp_conf.fai = fai_load(optarg);
              if (mplp_conf.fai == 0)  {
                   free(mplp_conf.fa);
                   return 1;
              }
              break;

         case 'c': 
              snvcall_conf.flag |= SNVCALL_CONS_AS_REF;
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
              mplp_conf.min_bq = atoi(optarg); 
              break;

         case 'Q': 
              snvcall_conf.min_altbq = atoi(optarg); 
              break;

         case 'a':
              snvcall_conf.def_altbq = atoi(optarg); 
              break;
/*
         case 'B': 
              mplp_conf.flag &= ~MPLP_REALN; 
              break;
*/
         case 'E': 
              mplp_conf.flag |= MPLP_REALN; /* BAQ */
              mplp_conf.flag |= MPLP_REDO_BAQ; /* ext BAQ */
              break;
              
         case 'm': 
              mplp_conf.min_mq = atoi(optarg); 
              break;

         case 'M': 
              mplp_conf.max_mq = atoi(optarg); 
              break;

         case 'J': 
              snvcall_conf.flag &= ~SNVCALL_USE_MQ; 
              break;

#ifdef USE_SOURCEQUAL
         case 'S': 
              snvcall_conf.flag &= ~SNVCALL_USE_SQ;
              break;
#endif

         case 's': 
              snvcall_conf.sig = strtof(optarg, (char **)NULL); /* atof */
              if (0==snvcall_conf.sig) {
                   LOG_FATAL("%s\n", "Couldn't parse sign-threshold"); 
                   return 1;
              }
              break;

         case 'b': 
              if (0 == strncmp(optarg, "auto", 4)) {
                   bonf_auto = 1;
                   snvcall_conf.bonf_dynamic = 0;

              } else if (0 == strncmp(optarg, "dynamic", 7)) {
                   snvcall_conf.bonf_dynamic = 1;
                   bonf_auto = 0;

              } else {
                   bonf_auto = 0;
                   snvcall_conf.bonf_dynamic = 0;

                   snvcall_conf.bonf = strtoll(optarg, (char **)NULL, 10); /* atol */ 
                   if (1>snvcall_conf.bonf) {
                        LOG_FATAL("%s\n", "Couldn't parse Bonferroni factor"); 
                        return 1;
                   }
              }
              break;

         case 'C': 
              snvcall_conf.min_cov = atoi(optarg); 
              break;
/*
         case 'd': 
              mplp_conf.max_depth = atoi(optarg); 
              break;
*/       
         case 'I': 
              mplp_conf.flag |= MPLP_ILLUMINA13; 
              break;

         /* already set: use-orphan, plp-summary-only, verbose, debug */

         case 'h': 
              usage(& mplp_conf, & snvcall_conf); 
              return 0; /* WARN: not printing defaults if some args where parsed */

         case '?': 
              LOG_FATAL("%s\n", "unrecognized arguments found. Exiting...\n"); 
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

    if (use_orphan) {
         mplp_conf.flag &= ~MPLP_NO_ORPHAN;
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


    /* FIXME: implement check_user_arg_logic() */
    assert(mplp_conf.min_mq <= mplp_conf.max_mq);
    assert(mplp_conf.min_bq <= snvcall_conf.min_altbq);
    assert(! (mplp_conf.bed && mplp_conf.reg));

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
              snvcall_conf.out = stdout;
         } else {
              if (NULL == (snvcall_conf.out = fopen(vcf_out, "w"))) {
                   LOG_FATAL("Couldn't open file '%s'. Exiting...\n", vcf_out);
                   return 1;
              }
         }
    } else {
         vcf_tmp_out = strdup(mktemp(vcf_tmp_template));
         if (NULL == vcf_tmp_out) {
              LOG_FATAL("%s\n", "Couldn't create temporary vcf file");
              return 1;
         }
         if (NULL == (snvcall_conf.out = fopen(vcf_tmp_out, "w"))) {
              LOG_FATAL("Couldn't open file '%s'. Exiting...\n", vcf_tmp_out);
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
    }


    if (bonf_auto && ! plp_summary_only) {
         if (! bed_file) {
              LOG_FATAL("%s\n", "Need bed-file for auto bonferroni correction.");
              return 1;
         }

#if 0
         double cov_mean;
         long long int num_non0cov_pos;
         LOG_DEBUG("Automatically determining Bonferroni factor for bam=%s reg=%s bed=%s\n",
                   bam_file, mplp_conf.reg, bed_file); 
         if (depth_stats(&cov_mean, &num_non0cov_pos, bam_file, mplp_conf.reg, bed_file,
                         &mplp_conf.minbq, &mplp_conf.min_mq)) {
              LOG_FATAL("%s\n", "Couldn't determine Bonferroni factor automatically\n"); 
              return 1;
         }
         snvcall_conf.bonf = num_non0cov_pos*3;
#else

         snvcall_conf.bonf = bonf_from_bedfile(bed_file);
         if (snvcall_conf.bonf<1) {
              LOG_FATAL("Automatically determining Bonferroni from bed"
                        " regions listed in %s failed\n", bed_file);
              return 1;
         }

#endif
         LOG_VERBOSE("Automatically determined Bonferroni factor = %lld\n", snvcall_conf.bonf);
    }

    if (debug) {
         dump_mplp_conf(& mplp_conf, stderr);
         dump_snvcall_conf(& snvcall_conf, stderr);
    }

    if (! plp_summary_only) {
         /* or use PACKAGE_STRING */
         vcf_write_header(snvcall_conf.out, mplp_conf.cmdline, mplp_conf.fa);
         plp_proc_func = &call_snvs;
    } else {
         plp_proc_func = &plp_summary;

    }

    rc = mpileup(&mplp_conf, plp_proc_func, (void*)&snvcall_conf,
                 1, (const char **) argv + optind + 1);


    if (snvcall_conf.out != stdout) {
         fclose(snvcall_conf.out);
    }

    /* snv calling completed. now filter according to the following schema:
     *  1. no_default_filter and ! dyn
     *     print
     *  2.1 no_ default_filter and dyn
     *     filter snvphred only
     *  2.2 ! no_default_filter and dyn
     *     filter snvphred and default
     *  2.3 ! no_default_filter and ! dyn 
     *     filter default
     */
    if (no_default_filter && ! snvcall_conf.bonf_dynamic) {
         /* vcf file needs no filtering and was already printed to
          * final destination. already taken care of above. */
         LOG_VERBOSE("%s\n", "No filtering needed. Variants already written to final destination");

    } else {
         char base_cmd[BUF_SIZE];
         char full_cmd[BUF_SIZE];
         snprintf(base_cmd, BUF_SIZE, 
                  "lofreq2_filter.py -p -i %s -o %s",
                  vcf_tmp_out, NULL==vcf_out ? "-" : vcf_out);

         if (no_default_filter && snvcall_conf.bonf_dynamic) {
              snprintf(full_cmd, BUF_SIZE, 
                      "%s --no-defaults --snv-phred %d", 
                      base_cmd, PROB_TO_PHREDQUAL(snvcall_conf.sig/snvcall_conf.bonf));

         } else if (! no_default_filter && snvcall_conf.bonf_dynamic) {
              snprintf(full_cmd, BUF_SIZE, 
                      "%s --snv-phred %d", 
                      base_cmd, PROB_TO_PHREDQUAL(snvcall_conf.sig/snvcall_conf.bonf));

         } else if (! no_default_filter && ! snvcall_conf.bonf_dynamic) {
              snprintf(full_cmd, BUF_SIZE, "%s", base_cmd);

         } else {
              LOG_FATAL("%s\n", "internal logic error during filtering");
              return 1;
         }

         LOG_VERBOSE("Executing %s\n", full_cmd);
         if (0 != (rc = system(full_cmd))) {
              LOG_ERROR("The following command failed: %s\n", full_cmd);
         } else {
              (void) unlink(vcf_tmp_out);
         }
    }

    free(vcf_tmp_out);
    free(vcf_out);
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
