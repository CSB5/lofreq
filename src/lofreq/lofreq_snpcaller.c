/* -*- c-file-style: "k&r" -*-
 *
 * This file is partially based on samtools' bam_plcmd.c
 *
 * FIXME missing license
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

#include "faidx.h"
#include "sam.h"
#include "kstring.h"

/* from bedidx.c */
void *bed_read(const char *fn);
void bed_destroy(void *_h);
int bed_overlap(const void *_h, const char *chr, int beg, int end);

#include "snpcaller.h"
#include "vcf.h"
#include "bam2depth.h"
#include "fet.h"
#include "utils.h"
#include "log.h"
#include "plp.h"
#include "bed.h"


#define SNVCALL_USE_MQ      0x10
#define SNVCALL_USE_SQ      0x20
#define SNVCALL_CONS_AS_REF 0x40

#define BUF_SIZE 1024


typedef struct {
     int min_altbq, def_altbq;/* tag:snvcall */
     long long int bonf;/* tag: snvcall */
     float sig;/* tag: snvcall */
     FILE *out;
     int flag;
} snvcall_conf_t;



void
report_var(FILE *stream, const plp_col_t *p, 
                const char ref, const char alt, 
                const float af, const int qual,
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



/* "Merge" MQ and BQ if requested and if MAQP not 255 (not available):
 *  P_jq = P_mq * + (1-P_mq) P_bq.
 */
int
merge_baseq_and_mapq(const int bq, const int mq)
{
     double mp, bp, jp; /* corresponding probs */
     int jq;

     if (mq == 255) {
          return bq;
     }
      
     /* No need to do computation in phred-space as
      * numbers won't get small enough.
      */
     mp = PHREDQUAL_TO_PROB(mq);
     bp = PHREDQUAL_TO_PROB(bq);

     jp = mp + (1.0 - mp) * bp;
     jq = PROB_TO_PHREDQUAL(jp);
#ifdef DEBUG
     LOG_DEBUG("P_M + (1-P_M) P_B:   %g + (1.0 - %g) * %g = %g  ==  Q%d + (1.0 - Q%d) * Q%d  =  Q%d\n",
               mp, mp, bp, jp, mq+33, mq+33, bq+33, jq+33);
#endif
     return jq;
}
/* merge_baseq_and_mapq() */



void
plp_summary(const plp_col_t *plp_col, const void* confp) 
{
     FILE* stream = stdout;
     int i;

     fprintf(stream, "%s\t%d\t%c\t%c",
             plp_col->target, plp_col->pos+1, plp_col->ref_base, plp_col->cons_base);
     for (i=0; i<NUM_NT4; i++) {
          fprintf(stream, "\t%c:%lu/%lu",
                  bam_nt4_rev_table[i],
                  plp_col->fw_counts[i],
                  plp_col->rv_counts[i]);
     }

     fprintf(stream, "\theads:%d\ttails:%d", plp_col->num_heads, plp_col->num_tails);
     fprintf(stream, "\tins=%d\tdels=%d", plp_col->num_ins, plp_col->num_dels);
     fprintf(stream, "\n");

     LOG_FIXME("%s\n", "unfinished");
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
call_lowfreq_snps(const plp_col_t *p, const void *confp)
{
     int *quals; /* qualities passed down to snpcaller */
     int quals_len; /* #elements in quals */
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
     if (p->num_dels || p->num_ins) {
          LOG_FIXME("%s:%d (p->num_dels=%d p->del_quals=%d p->num_ins=%d p->ins_quals.n=%d\n", 
                    p->target, p->pos+1, p->num_dels, p->del_quals.n, p->num_ins, p->ins_quals.n);
          if (p->num_dels && p->del_quals.n) {
               LOG_FIXME("Call deletions at %s:%d\n", p->target, p->pos+1);
          }
          if (p->num_ins && p->ins_quals.n) {
               LOG_FIXME("Call insertions at %s:%d\n", p->target, p->pos+1);
          }
     }

     /* CONSVAR, i.e. the consensus determined here is different from
      * the reference coming from a fasta file 
      */
     if (p->ref_base != p->cons_base && !cons_as_ref && p->ref_base != 'N') {
          const int is_indel = 0;
          const int is_consvar = 1;
          const int qual = -1;
          float af = (p->fw_counts[bam_nt4_table[(int) p->cons_base]] 
                      + 
                      p->rv_counts[bam_nt4_table[(int) p->cons_base]]) 
               / (float)p->coverage;

          report_var(conf->out, p, p->ref_base, p->cons_base, af, qual, is_indel, is_consvar);
          LOG_DEBUG("cons var snp: %s %d %c>%c\n",
                    p->target, p->pos+1, p->ref_base, p->cons_base);          
     }

     if (NULL == (quals = malloc(p->coverage * sizeof(int)))) {
          /* coverage = base-count after read level filtering */
          fprintf(stderr, "FATAL: couldn't allocate memory at %s:%s():%d\n",
                  __FILE__, __FUNCTION__, __LINE__);
          free(quals);
          return;
     }
    
     quals_len = 0;
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
               int bq, mq, final_q;
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
                    final_q = merge_baseq_and_mapq(bq, mq);

               } else {
                    final_q = bq;
               }

               quals[quals_len++] = final_q;
          }
     }

     for (i=0; i<3; i++) {
          if (alt_counts[i]) {
               got_alt_bases = 1;
               break;
          }
     }
     if (! got_alt_bases) {
          LOG_DEBUG("%s %d: only cons bases left after filtering.\n", p->target, p->pos+1);
          /* ...and CONSVAR already reported */
          free(quals);
          return;
     }

     /* sorting in theory should be numerically more stable and also
      * make snpcallerfaster */
     qsort(quals, quals_len, sizeof(int), int_cmp);

     LOG_DEBUG("%s %d: passing down %d quals with noncons_counts (%d, %d, %d) to snpcaller()\n",
               p->target, p->pos+1, quals_len, alt_counts[0], alt_counts[1], alt_counts[2]);

     if (snpcaller(pvalues, quals, quals_len, 
                  alt_counts, conf->bonf, conf->sig)) {
          fprintf(stderr, "FATAL: snpcaller() failed at %s:%s():%d\n",
                  __FILE__, __FUNCTION__, __LINE__);
          free(quals);
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
               continue;
          }
          if (p->ref_base==alt_base && !cons_as_ref) {
               continue;
          }

          if (pvalue * (double)conf->bonf < conf->sig) {
               const int is_indel = 0;
               const int is_consvar = 0;
               float af = alt_raw_count/(float)p->coverage;

               report_var(conf->out, p, reported_snv_ref, alt_base, 
                          af, PROB_TO_PHREDQUAL(pvalue), 
                          is_indel, is_consvar);
               LOG_DEBUG("low freq snp: %s %d %c>%c pv-prob:%g;pv-qual:%d counts-raw:%d/%d=%.6f counts-filt:%d/%d=%.6f\n",
                         p->target, p->pos+1, p->cons_base, alt_base,
                         pvalue, PROB_TO_PHREDQUAL(pvalue),
                         /* counts-raw */ alt_raw_count, p->coverage, alt_raw_count/(float)p->coverage,
                         /* counts-filt */ alt_count, quals_len, alt_count/(float)quals_len);
          }
     }
     free(quals);
}
/* call_lowfreq_snps() */


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

     LOG_FIXME("src_pvalue = %g. Actual prob = %g\n", src_pvalue, exp(probvec[n_mismatches]));

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



void
dump_snvcall_conf(const snvcall_conf_t *c, FILE *stream) 
{
     fprintf(stream, "snvcall options\n");
     fprintf(stream, "  min_altbq = %d\n", c->min_altbq);
     fprintf(stream, "  def_altbq = %d\n", c->def_altbq);
     fprintf(stream, "  bonf         = %lld  (might get recalculated later)\n", c->bonf);
     fprintf(stream, "  sig          = %f\n", c->sig);
     fprintf(stream, "  out          = %p\n", (void*)c->out);
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
     /* generic */
     fprintf(stderr, "           --verbose                be verbose\n");
     fprintf(stderr, "           --debug                  enable debugging\n");
     /* regions */                                        
     fprintf(stderr, "       -r | --region STR            region in which pileup should be generated [null]\n");
     fprintf(stderr, "       -l | --bed FILE              list of positions (chr pos) or regions (BED) [null]\n");
     /*  */                                               
     fprintf(stderr, "       -d | --maxdepth INT          max per-BAM depth to avoid excessive memory usage [%d]\n", mplp_conf->max_depth);
     fprintf(stderr, "       -f | --reffa FILE            faidx indexed reference sequence file [null]\n");
     fprintf(stderr, "       -c | --cons-as-ref           Use consensus base as ref, i.e. ignore base given in reffa (reffa still used for BAQ)\n");
     /* */                                                
     fprintf(stderr, "       -o | --out FILE              vcf output file [- = stdout]\n");
     /* base call quality and baq */                      
     fprintf(stderr, "       -q | --min-bq INT            skip any base with baseQ smaller than INT [%d]\n", mplp_conf->min_bq);
     fprintf(stderr, "       -Q | --min-altbq INT         skip nonref-bases with baseQ smaller than INT [%d]. Not active if ref is N\n", snvcall_conf->min_altbq);
     fprintf(stderr, "       -a | --def-altbq INT         nonref base qualities will be replace with this value [%d]\n", snvcall_conf->def_altbq);
     fprintf(stderr, "       -B | --no-baq                disable BAQ computation\n");
     /* fprintf(stderr, "       -E           extended     BAQ for higher sensitivity but lower specificity\n"); */
     /* mapping quality */                                
     fprintf(stderr, "       -m | --min_mq INT            skip alignments with mapQ smaller than INT [%d]\n", mplp_conf->min_mq);
     fprintf(stderr, "       -M | --max_mq INT            cap mapping quality at INT [%d]\n", mplp_conf->max_mq);
     fprintf(stderr, "       -J | --no-mq                 don't merge mapQ into baseQ: P_e = P_mq + (1-P_mq) P_bq\n");
#ifdef USE_SOURCEQUAL                                     
     fprintf(stderr, "       -S | --no-sq                 don't merge sourceQ into baseQ\n");
#endif                                                    
     /* stats */                                          
     fprintf(stderr, "       -s | --sig                   P-value cutoff / significance level [%f]\n", snvcall_conf->sig);
     fprintf(stderr, "       -b | --bonf                  Bonferroni factor. INT or 'auto' (needs bed-file) [auto]\n");
     /* misc */                                           
     fprintf(stderr, "       -I | --illumina-1.3          assume the quality is Illumina-1.3-1.7/ASCII+64 encoded\n");
     fprintf(stderr, "            --use-orphan            count anomalous read pairs\n");
     fprintf(stderr, "            --plp-summary-only      no snv-calling. just output pileup summary per column\n");
     fprintf(stderr, "       -p | --pseudo-parallel INT   pseudo-parallel using INT processors, by splitting up the bed-file. NOTE: don't use if disk IO is a bottleneck\n");
}
/* usage() */


/* FIXME evil hack, use MP/I in main routine instead */
int 
main_call_pseudo_parallel(int argc, char *argv[], const int num_proc, 
                          const char *bed_file,  FILE *vcfout)
{
     int i;
     long int line_count;
     int num_splits;
     char *lofreq_args = NULL;
     char *tmpdir;
     /* FIXME hardcoded /tmp/ */
     char templ[] = "/tmp/lofreq2-call-pseudoparallel.XXXXXX";
     char cmd_buf[BUF_SIZE];
     const char out_xargs[] = "-o @.vcf \0";
     const char bed_xargs[] = "-l @ \0";
     char **split_vcf_files;
     int num_split_vcf_files;

     if (! bed_file) {
          LOG_FATAL("%s\n", "Pseudo-parallel mode requires a bed-file");
          return -1;
     }
     if (num_proc<2) {
          LOG_FATAL("%s\n", "Running in pseudo-parallel mode with less"
                    " then 2 processors doesn't make sense");
          return -1;
     }

     /* FIXME this will also count empty lines */
     line_count =  count_lines(bed_file);
     if (line_count<2) {
          LOG_ERROR("Your bed-file has only %d entries."
                    " Can't run in pseudo-parallel mode.\n", 
                    line_count);
          return -1;        
     }
         
     LOG_VERBOSE("Will try to run %d split processes\n", num_proc);

     /* split in how many files */
     if (line_count<num_proc) {
          LOG_WARN("The number of regions in your bed-file is smaller than"
                   " the number of requested pseudo-parallel threads."
                   " That's fine, but I can't use all %d threads\n", num_proc);             
          num_splits = line_count;
     } else {
          num_splits = line_count / num_proc;
     }

     /* use system command 'split' to split files */
     if (NULL == (tmpdir = mkdtemp(templ))) {
          LOG_FATAL("%s\n", "Couldn't create temporary directory");
          return -1;
     }
     if (0 > snprintf(cmd_buf, BUF_SIZE, "split -l %d %s %s/",
                      num_splits, bed_file, tmpdir)) {
          LOG_FATAL("%s\n", "snprintf failed");
          return -1;
     }
     LOG_DEBUG("Running: %s\n", cmd_buf);
     if (system(cmd_buf)) {
          LOG_FATAL("%s\n", "Splitting your bed-file with split failed.");
          return -1;
     }
         
     lofreq_args = calloc(1, sizeof(char));
     for (i=0; i<argc-1 /* leave out BAM */; i++) {
          if (0 == strcmp(argv[i], "--pseudo-parallel") || 0 == strcmp(argv[i], "-p")) {
               i+=1; /* skip arg */ 
               continue;
          }
          if (0 == strcmp(argv[i], "--out") || 0 == strcmp(argv[i], "-o")) {
               i+=1; /* skip arg */ 
               continue;
          }
          if (0 == strcmp(argv[i], "--bed") || 0 == strcmp(argv[i], "-l")) {
               i+=1; /* skip arg */ 
               continue;
          }
          lofreq_args = realloc(lofreq_args,
                                strlen(lofreq_args) + strlen(argv[i]) + 1/*space*/); /* FIXME inefficient */
          lofreq_args = strcat(lofreq_args, argv[i]);
          lofreq_args = strcat(lofreq_args, " ");
     }
     lofreq_args = realloc(lofreq_args,
                           strlen(lofreq_args) + strlen(out_xargs) + 1);
     lofreq_args = strcat(lofreq_args, out_xargs);
     
     lofreq_args = realloc(lofreq_args,
                           strlen(lofreq_args) + strlen(bed_xargs) + 1);
     lofreq_args = strcat(lofreq_args, bed_xargs);
     
     /* BAM */
     lofreq_args = realloc(lofreq_args,
                           strlen(lofreq_args) + strlen(argv[argc-1]) + 1/*space*/); /* FIXME inefficient */
     lofreq_args = strcat(lofreq_args, argv[argc-1]);
     
     /* FIXME: replace xargs with forks */

     /* asterisk needed, otherwise @ becomes just basename */
     snprintf(cmd_buf, BUF_SIZE, "ls %s/* | grep -v vcf$ | xargs -I@ -n 1 -P %d %s",
              tmpdir, num_proc, lofreq_args);
     
     LOG_FIXME("Running: %s\n", cmd_buf);
     /* FIXME catch err by redirecting stderr to tmpfile */
     /* could try popen her to parse output directly, but all
      * subprocesses will print a header which might get mingled.
      */
     if (system(cmd_buf)) {
          LOG_FATAL("%s\n", "Running in pseudo-parallel mode failed (at xargs stage)");
          return -1;
     }
     

     num_split_vcf_files = ls_dir(&split_vcf_files, tmpdir, ".vcf", 1);
     if (num_split_vcf_files < 1) {
          LOG_FATAL("%s\n", "Couldn't list output vcf files");
          return -1;
     }
         
     /* the sub vcf's are sorted lexicographically, which also means
      * in bed order. therefore only print header from first.
      */
     for (i=0; i<num_split_vcf_files; i++) {
          char *fn = split_vcf_files[i];
          FILE *fh;
          char line[BUF_SIZE];
          LOG_FIXME("Reading from %s\n", fn);
          if (NULL == (fh = fopen(fn, "r"))) {
               LOG_FATAL("Couldn't open vcf file %s", fn);
               return -1;
          }
          while (NULL != fgets(line, BUF_SIZE, fh)) {
               if ((0==i && line[0]=='#') || (line[0]!='#')) {
                    fprintf(vcfout, "%s", line);
               }
          }
          fclose(fh);
     }

     for (i=0; i<num_split_vcf_files; i++) {
          free(split_vcf_files[i]);
     }
     free(split_vcf_files);
     free(lofreq_args);
     if (0) {
          (void) unlink(tmpdir);
     } else {
          LOG_FIXME("Not unlinking tmpdir %s\n", tmpdir);
     }
     return 0;
}


int 
main_call(int argc, char *argv[])
{
     /* based on bam_mpileup() */
     int c, i;
     static int use_orphan = 0;
     static int plp_summary_only = 0;
     int bonf_auto = 1;
     char *bam_file;
     char *bed_file = NULL;
     bed_t bed;
     mplp_conf_t mplp_conf;
     snvcall_conf_t snvcall_conf;
     /*void (*plp_proc_func)(const plp_col_t*, const snvcall_conf_t*);*/
     void (*plp_proc_func)(const plp_col_t*, const void*);
     int pseudo_parallel = 0;
     int rc = 0;

     for (i=0; i<argc; i++) {
          LOG_DEBUG("arg %d: %s\n", i, argv[i]);
     }
     LOG_FIXME("%s\n", "- Proper source qual use missing");
     LOG_FIXME("%s\n", "- Indel handling missing");
     LOG_FIXME("%s\n", "- Implement routine test against old SNV caller");
     LOG_FIXME("%s\n", "- Test actual SNV and SB values for both types of SNVs");


     memset(&bed, 0, sizeof(bed_t));

     /* default pileup options */
     memset(&mplp_conf, 0, sizeof(mplp_conf_t));
     mplp_conf.max_mq = 255;
     mplp_conf.min_bq = 3;
     mplp_conf.capQ_thres = 0;
     mplp_conf.max_depth = 1000000;
     mplp_conf.flag = MPLP_NO_ORPHAN | MPLP_REALN | MPLP_EXT_BAQ;
    
     /* default snvcall options */
     memset(&snvcall_conf, 0, sizeof(snvcall_conf_t));
     snvcall_conf.min_altbq = 20;
     snvcall_conf.def_altbq = snvcall_conf.min_altbq;
     snvcall_conf.bonf = 1;
     snvcall_conf.sig = 0.05;
     snvcall_conf.out = stdout;
     snvcall_conf.flag = SNVCALL_USE_MQ;/* | MPLP_USE_SQ; FIXME */

    /* keep in sync with long_opts_str and usage */
    while (1) {
         static struct option long_opts[] = {
              /* see usage sync */
              {"verbose", no_argument, &verbose, 1},
              {"debug", no_argument, &debug, 1},

              {"region", required_argument, NULL, 'r'},
              {"bed", required_argument, NULL, 'l'}, /* changes here must be reflected in pseudo_parallel code as well */
              
              {"maxdepth", required_argument, NULL, 'd'},
              {"reffa", required_argument, NULL, 'f'},
              {"cons-is-ref", no_argument, NULL, 'c'},

              {"out", required_argument, NULL, 'o'}, /* NOTE changes here must be reflected in pseudo_parallel code as well */

              {"min-bq", required_argument, NULL, 'q'},
              {"min-altbq", required_argument, NULL, 'Q'},
              {"def-altbq", required_argument, NULL, 'a'},
              {"no-baq", no_argument, NULL, 'B'},
              /*{"ext-baq", required_argument, NULL, 'E'},*/
                   
              {"min-mq", required_argument, NULL, 'm'},
              {"max-mq", required_argument, NULL, 'M'},
              {"no-mq", no_argument, NULL, 'J'},
#ifdef USE_SOURCEQUAL
              {"no-sq", no_argument, NULL, 'S'},
#endif
              {"bonf", required_argument, NULL, 'b'},
              {"sig", required_argument, NULL, 's'},
                   
              {"illumina-1.3", no_argument, NULL, 'I'},
              {"use-orphan", no_argument, &use_orphan, 1},
              {"plp-summary-only", no_argument, &plp_summary_only, 1},
              {"pseudo-parallel", required_argument, NULL, 'p'}, /* name change here has to be reflected further down in code as well */

              {"help", no_argument, NULL, 'h'},

              {0, 0, 0, 0} /* sentinel */
         };

         /* keep in sync with long_opts and usage */
         static const char *long_opts_str = "r:l:d:f:co:q:Q:a:Bm:M:JSb:s:Iph"; 

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

         case 'c': 
              snvcall_conf.flag |= SNVCALL_CONS_AS_REF;
              break;

         case 'l': 
              /*mplp_conf.bed = bed_read(optarg); */
              bed_file = strdup(optarg);
              break;

         case 'd': 
              mplp_conf.max_depth = atoi(optarg); 
              break;

         case 'f':
              mplp_conf.fa = strdup(optarg);
              mplp_conf.fai = fai_load(optarg);
              if (mplp_conf.fai == 0) 
                   return 1;
              break;

         case 'o':
              if (0 != strcmp(optarg, "-")) {
                   if (file_exists(optarg)) {
                        LOG_FATAL("Cowardly refusing to overwrite file '%s'. Exiting...\n", optarg);
                        return 1;
                   }
                   if (NULL == (snvcall_conf.out = fopen(optarg, "w"))) {
                        LOG_FATAL("Couldn't open file '%s'. Exiting...\n", optarg);
                        return 1;
                   }
              } /* else: already set to stdout */
              break;

         case 'q': 
              LOG_FIXME("optarg=%s\n", optarg);
              mplp_conf.min_bq = atoi(optarg); 
              break;

         case 'Q': 
              snvcall_conf.min_altbq = atoi(optarg); 
              break;

         case 'a': snvcall_conf.def_altbq = atoi(optarg); 
              break;

         case 'B': 
              mplp_conf.flag &= ~MPLP_REALN; 
              break;
         /* case 'E': mplp.flag |= MPLP_EXT_BAQ; break; */
              
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
         case 'I': 
              mplp_conf.flag |= MPLP_ILLUMINA13; 
              break;

         case 'b': 
              if (0 == strncmp(optarg, "auto", 4)) {
                   bonf_auto = 1;

              } else {
                   bonf_auto = 0;
                   snvcall_conf.bonf = strtoll(optarg, (char **)NULL, 10); /* atol */ 
                   if (1>snvcall_conf.bonf) {
                        LOG_FATAL("%s\n", "Couldn't parse Bonferroni factor"); 
                        return 1;
                   }
              }
              break;

         case 's': 
              snvcall_conf.sig = strtof(optarg, (char **)NULL); /* atof */
              if (0==snvcall_conf.sig) {
                   LOG_FATAL("%s\n", "Couldn't parse sign-threshold"); 
                   return 1;
              }
              break;
              
         case 'p':
              /* FIXME if (NULL==optarg) {argh};e.g. -(!)pseudparallel */
              pseudo_parallel = atoi(optarg); 
              break;

         case 'h': 
              usage(& mplp_conf, & snvcall_conf); 
              return 0; /* WARN: not printing defaults if some args where parsed */

         case '?': 
              LOG_FATAL("%s\n", "unrecognized arguments found. Exiting...\n"); 
              return 1;
#if 0
         case 0:
              fprintf(stderr, "ERROR: long opt (%s) not mapping to short option. Exiting...\n", long_opts[long_opts_index].name); 
              return 1;
#endif
         default:
              break;
         }
    }

    if (use_orphan) {
         mplp_conf.flag &= ~MPLP_NO_ORPHAN;
    }

    if ( ! (snvcall_conf.flag & SNVCALL_CONS_AS_REF) && ! mplp_conf.fa && ! plp_summary_only) {
         LOG_FATAL("%s\n", "Need a reference when not calling in consensus mode...\n"); 
         return 1;
    }
    if (! plp_summary_only & ! mplp_conf.fa) {
         LOG_WARN("%s\n", "Calling SNVs without reference\n"); 
    }

    mplp_conf.cmdline[0] = '\0';
    for (i=0; i<argc; i++) {
         strncat(mplp_conf.cmdline, argv[i], sizeof(mplp_conf.cmdline)-strlen(mplp_conf.cmdline)-2);
         strcat(mplp_conf.cmdline, " ");
    }

    if (argc == 2) {
        fprintf(stderr, "\n");
        usage(& mplp_conf, & snvcall_conf);
        return 1;
    }
    if (1 != argc - optind - 1) {
        fprintf(stderr, "Need exactly one BAM file as last argument\n");
        return 1;
    }
    bam_file = (argv + optind + 1)[0];
    if (0 == strcmp(bam_file, "-") && mplp_conf.reg) {
         LOG_FATAL("%s\n", "Need index if region was given and index file can't be provided when using stdin mode.");
         return 1;         
    }

    if (bed_file) {
         /* FIXME should be using bed_read and mplp.bed instead of
          * parsing bed file again */
         if (-1 == parse_bed(&bed, bed_file)) {
              LOG_FATAL("Parsing of %s failed\n", bed_file); 
              return 1;
         }
#if 1
         LOG_FIXME("%s\n", "debug bed dumping"); 
         dump_bed(&bed);
#endif
    }

    if (bonf_auto && ! plp_summary_only) {
         double cov_mean;
         long long int num_non0cov_pos;

         if (pseudo_parallel>1) {
              LOG_FATAL("%s\n", "Can't use auto bonf in pseudo-parallel mode.");
              return 1;
         }
         if (! bed_file) {
              LOG_FATAL("%s\n", "Need bed-file for auto bonferroni correction.");
              return 1;
         }
         if (plp_summary_only) {
              LOG_FATAL("%s\n", "Automatically determining a Bonferroni"
                        " factor for printing a pileup summary only"
                        " doesn't make sense.");
              return 1;
         }

#if 0
         LOG_DEBUG("Automatically determining Bonferroni factor for bam=%s reg=%s bed=%s\n",
                   bam_file, mplp_conf.reg, bed_file); 
         if (depth_stats(&cov_mean, &num_non0cov_pos, bam_file, mplp_conf.reg, bed_file,
                         &mplp_conf.min_bq, &mplp_conf.min_mq)) {
              LOG_FATAL("%s\n", "Couldn't determine Bonferroni factor automatically\n"); 
              return 1;
         }
         snvcall_conf.bonf = num_non0cov_pos*3;
#else

         snvcall_conf.bonf = 3*bed_pos_sum(&bed);
         if (snvcall_conf.bonf<1) {
              LOG_FATAL("Automatically determining Bonferroni from bed"
                        " regions listed in %s failed\n", bed_file);
         }

#endif
         LOG_VERBOSE("Automatically determined Bonferroni factor = %lld\n", snvcall_conf.bonf);
    }

    if (debug) {
         dump_mplp_conf(& mplp_conf, stderr);
         dump_snvcall_conf(& snvcall_conf, stderr);
    }

    /* FIXME: implement logic_check_opts() */
    assert(mplp_conf.min_mq <= mplp_conf.max_mq);
    assert(mplp_conf.min_bq <= snvcall_conf.min_altbq);
    assert(! (mplp_conf.bed && mplp_conf.reg));

    if (pseudo_parallel>1) { /* to be able to break, which wouldn't be possible with an if */
         rc = main_call_pseudo_parallel(argc, argv, pseudo_parallel, 
                                        bed_file, snvcall_conf.out);
         goto cleanup; 
    }
   
    if (! plp_summary_only) {
         /* FIXME would be nice to use full command line here instead of PACKAGE_STRING */
         vcf_write_header(snvcall_conf.out, PACKAGE_STRING, mplp_conf.fa);
         plp_proc_func = &call_lowfreq_snps;
    } else {
         plp_proc_func = &plp_summary;

    }

    if (bed.nregions) {
         for (i=0; i<bed.nregions; i++) {
              char buf[BUF_SIZE];
              snprintf(buf, BUF_SIZE, "%s:%d-%d", bed.region[i].chrom, bed.region[i].start, bed.region[i].end);
              mplp_conf.reg = strdup(buf);
              LOG_VERBOSE("Running on region %s with bonf %lld (mplp_conf.bed is %p)\n", 
                          mplp_conf.reg, snvcall_conf.bonf, mplp_conf.bed);
              (void) mpileup(&mplp_conf, plp_proc_func, (void*)&snvcall_conf,
                             1, (const char **) argv + optind + 1);
              free(mplp_conf.reg);
              mplp_conf.reg = NULL;
         }
    } else {
         /* no bed. no reg. run on the whole BAM file */
         (void) mpileup(&mplp_conf, plp_proc_func, (void*)&snvcall_conf,
                        1, (const char **) argv + optind + 1);
    }

    rc = 0;

cleanup:

    if (snvcall_conf.out != stdout) {
         fclose(snvcall_conf.out);
    }
    LOG_FIXME("%s\n", "before mplp_conf free");
    free(mplp_conf.reg); 
    free(mplp_conf.fa);
    LOG_FIXME("%s\n", "before fai destroy");
    if (mplp_conf.fai) {
         fai_destroy(mplp_conf.fai);
    }
    LOG_FIXME("%s\n", "before free_bed");
    free_bed(&bed);
    free(bed_file);
/*
    if (mplp_conf.bed) {
         bed_destroy(mplp_conf.bed);
    }
*/
    LOG_VERBOSE("%s\n", "Successful exit.");
    return rc;
}
/* main_call */


#ifdef MAIN_TEST


int main()
{     
     fprintf(stderr, "WARNING: main replaced with test function. just monkeying around\n");

     if (1) {
          char test_nucs[] = "ACGTNRYacgtnryZ\0";
          int i;
          
          for (i=0; i<strlen(test_nucs); i++) {
               printf("%d %c - %d - %d - %d\n", i, test_nucs[i],
                      bam_nt16_table[(int)test_nucs[i]],
                      bam_nt16_nt4_table[bam_nt16_table[(int)test_nucs[i]]],
                      bam_nt4_table[(int)test_nucs[i]]);
          }
     }

     if (1) {
          int quals[] = {30, 30, 30, 30, 30, 30, 30, 30, 30, 30,
                         30, 30, 30, 30, 30, 30, 30, 30, 30, 30,
                         30, 30, 30, 30, 30, 30, 30, 30, 30, 30,
                         30, 30, 30, 30, 30, 30, 30, 30, 30, 30,
                         30, 30, 30, 30, 30, 30, 30, 30, 30, 30};
          int n_quals = 50;
          int n_mismatches = 1;
          double *probvec;
          double src_pvalue;
          int src_qual; 
          int i;

          qsort(quals, n_quals, sizeof(int), int_cmp);

          for (n_mismatches=0; n_mismatches<n_quals/2; n_mismatches++) {
               probvec = poissbin(&src_pvalue, quals,
                                  n_quals, n_mismatches, 1.0, 0.05);
               
               if (src_pvalue>1.0) {/* DBL_MAX is default return value */
                    src_pvalue = 1.0;/*-DBL_EPSILON;*/
               }
               /* src_pvalue: what's the chance of seeing n_mismatches or more
                * given quals? or: how likely is this read from the genome.
                * 1-src_value = prob read is not from genome
                */
               LOG_FIXME("Orig src_pv = %f", src_pvalue);
               src_pvalue = 1.0-src_pvalue;
               free(probvec);
               
               src_qual =  PROB_TO_PHREDQUAL(src_pvalue);
               fprintf(stderr, "| src_pv = %f = Q%d for %d/%d mismatches. All quals: ", 
                       src_pvalue, src_qual, n_mismatches, n_quals);
               for (i=0; i<n_quals; i++) {
                    fprintf(stderr, " %d", quals[i]);
               }
               fprintf(stderr, "\n");
          }
     }
     return 0;
}

#endif

