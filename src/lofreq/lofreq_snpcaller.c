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

#include "sam.h"
#include "faidx.h"
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


/* mpileup configuration flags 
 */
#define MPLP_NO_ORPHAN   0x10
#define MPLP_REALN       0x20
#define MPLP_EXT_BAQ     0x40
#define MPLP_ILLUMINA13  0x80
#define MPLP_IGNORE_RG   0x100
#define MPLP_USE_MQ      0x200
#define MPLP_USE_SQ      0x400


#define MYNAME "lofreq call"


const char *bam_nt4_rev_table = "ACGTN"; /* similar to bam_nt16_rev_table */

unsigned char bam_nt4_table[256] = {
     4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
     4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
     4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
     4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
     4,0,4,1,4,4,4,2,4,4,4,4,4,4,4,4,
     4,4,4,4,3,4,4,4,4,4,4,4,4,4,4,4,
     4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
     4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
     4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
     4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
     4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
     4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
     4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
     4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
     4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
     4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
};
#define NUM_NT4 5 /* strlen(bam_nt4_rev_table); */



/* mpileup configuration structure 
 * FIXME should be logically split into pileup, snvcall and others
 */
typedef struct {
     int max_mq, min_mq;
     int flag;
     int capQ_thres;
     int max_depth;
     int min_bq, min_altbq;
     int def_altbq;
     unsigned long int bonf;
     float sig;
     char *reg;
     char *fa;
     faidx_t *fai;
     void *bed;

     char cmdline[1024];
     FILE *out;
} mplp_conf_t;

typedef struct {
     bamFile fp;
     bam_iter_t iter;
     bam_header_t *h;
     int ref_id;
     char *ref;
     const mplp_conf_t *conf;
} mplp_aux_t;

typedef struct {
    int n;
    int *n_plp, *m_plp;
    bam_pileup1_t **plp;
} mplp_pileup_t;



typedef struct {
     char *target; /* chromsome or sequence name */
     int pos; /* position */
     char ref_base; /* uppercase reference base (given by fasta) */
     char cons_base; /* uppercase consensus base according to base-counts, after read-level filtering. */
     int coverage; /* coverage after read-level filtering i.e. same as in samtools mpileup (n_plp) but without indels! */

     /* list of qualities: keeping them all here in one place so that
      * filtering can become separate step. alternative is to filter
      * during pileup. the latter doesn't work if you want to filter
      * based on a consensus which you don't know in advance */
     int_varray_t base_quals[NUM_NT4]; 
     int_varray_t map_quals[NUM_NT4]; 
     int_varray_t source_quals[NUM_NT4]; 
     long int fw_counts[NUM_NT4]; 
     long int rv_counts[NUM_NT4]; 
     /* fw_counts[b] + rv_counts[b] = x_quals.n = coverage */

     int num_heads; /* number of read starts at this pos */
     int num_tails; /* number of read ends at this pos */

     int num_ins, sum_ins;
     int_varray_t ins_quals; 

     int num_dels, sum_dels;
     int_varray_t del_quals; 

     /* changes here should be reflected in plp_col_init, plp_col_free etc. */
} plp_col_t;


#define PLP_COL_ADD_QUAL(p, q)   int_varray_add_value((p), (q))


void plp_col_init(plp_col_t *p) {
    int i;

    p->target =  NULL;
    p->pos = -INT_MAX;
    p->ref_base = '\0';
    p->cons_base = 'N';
    p->coverage = -INT_MAX;
    for (i=0; i<NUM_NT4; i++) {
         int_varray_init(& p->base_quals[i], 0);
         int_varray_init(& p->map_quals[i], 0);
         int_varray_init(& p->source_quals[i], 0); /* FIXME unused */
         p->fw_counts[i] = 0;
         p->rv_counts[i] = 0;
    }

    p->num_heads = p->num_tails = 0;

    p->num_ins = p->sum_ins = 0;
    int_varray_init(& p->ins_quals, 0);
    p->num_dels = p->sum_dels = 0;
    int_varray_init(& p->del_quals, 0);
}

void plp_col_free(plp_col_t *p) {
    int i;

    free(p->target);
    for (i=0; i<NUM_NT4; i++) {
         int_varray_free(& p->base_quals[i]);
         int_varray_free(& p->map_quals[i]);
         int_varray_free(& p->source_quals[i]);
    }
    int_varray_init(& p->ins_quals, 0);
    int_varray_init(& p->del_quals, 0);
}

void plp_col_debug_print(const plp_col_t *p, FILE *stream)
{
     int i;
     
     fprintf(stream, "%s\t%d\t%c\t%c\tcounts:rv/fw",
             p->target, p->pos+1, p->ref_base, p->cons_base);
     for (i=0; i<NUM_NT4; i++) {
          fprintf(stream, " %c:%lu/%lu",
                  bam_nt4_rev_table[i],
                  p->fw_counts[i],
                  p->rv_counts[i]);
     }

     fprintf(stream, " heads:%d tails:%d", p->num_heads, p->num_tails);
     fprintf(stream, " ins=%d del=%d", p->num_ins, p->num_dels);
     fprintf(stream, "\n");

#if 0
     for (i=0; i<NUM_NT4; i++) {
          int j;
          fprintf(stream, "%c BQs (%lu): " , bam_nt4_rev_table[i], p->base_quals[i].n);
          for (j=0; j<p->base_quals[i].n; j++) {
               fprintf(stream, " %d", p->base_quals[i].data[j]);
          }
          fprintf(stream, "\n");
     }
#endif
}

/* attempt to keep a function in here that produces output similar to
 * the last pre-c version which can be easily parsed from Python. Note
 * however, that defaults have changed and that filtering was done differently before.
 */
void plp_col_mpileup_print(const plp_col_t *p, const mplp_conf_t *conf, FILE *stream)
{
     int i, j;
     
     fprintf(stream, "%s\t%d\t%c\t%d\t", 
             p->target, p->pos+1, p->ref_base,p->coverage);
     for (i=0; i<NUM_NT4; i++) {
          for (j=0; j<p->base_quals[i].n; j++) {
               fprintf(stream, "%c%c",
                       bam_nt4_rev_table[i],  p->base_quals[i].data[j]+33);
          }
     }
           
     fprintf(stream, "\t#heads=%d #tails=%d #ins=%d ins_len=%.1f #del=%d del_len=%.1f\n",
          p->num_heads, p->num_tails,
          p->num_ins, p->num_ins ? p->sum_ins/(float)p->num_ins : 0,
          p->num_dels, p->num_dels ? p->sum_dels/(float)p->num_dels : 0);
}



void report_var(FILE *stream, const plp_col_t *p, 
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


/* "Merge" MQ and BQ if requested and if MAQP not 255 (not available):
 *  P_jq = P_mq * + (1-P_mq) P_bq.
 */
int merge_baseq_and_mapq(const int bq, const int mq)
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


/* low-freq vars always called against cons_base, which might be
 * different from ref_base. if cons_base != ref_base then it's a
 * cons-var.
 * 
 * Assuming conf->min_bq and read-level filtering was already done
 * upstream. altbase mangling happens here however.
 * 
 */
void call_lowfreq_snps(const plp_col_t *p, const mplp_conf_t *conf)
{
     int *quals; /* qualities passed down to snpcaller */
     int quals_len; /* #elements in quals */
     int i, j;

     /* 4 bases ignoring N, -1 reference/consensus base makes 3 */
     double pvalues[3]; /* pvalues reported back from snpcaller */
     int alt_counts[3]; /* counts for alt bases handed down to snpcaller */
     int alt_raw_counts[3]; /* raw, unfiltered alt-counts */
     int alt_bases[3];/* actual alt bases */
     int alt_idx;
     int got_alt_bases = 0;

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

     /* check for consensus snps, i.e. those where the consensus
      * determined here is different from the reference coming from a
      * fasta file */
     if (p->ref_base != 'N' && p->ref_base != p->cons_base) {
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
               int bq, mq, sq, final_q;

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

               if ((conf->flag & MPLP_USE_MQ)) {
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


     for (i=0; i<3; i++) {
          int alt_base = alt_bases[i];
          int alt_count = alt_counts[i];
          int alt_raw_count = alt_raw_counts[i];
          double pvalue = pvalues[i];
          if (pvalue * (double)conf->bonf < conf->sig) {
               const int is_indel = 0;
               const int is_consvar = 0;
               float af = alt_raw_count/(float)p->coverage;
               report_var(conf->out, p, p->cons_base, alt_base, 
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


/* FIXME get rid in the future */
static inline int printw(int c, FILE *fp)
{
    char buf[16];
    int l, x;
    if (c == 0) return fputc('0', fp);
    for (l = 0, x = c < 0? -c : c; x > 0; x /= 10) buf[l++] = x%10 + '0';
    if (c < 0) buf[l++] = '-';
    buf[l] = 0;
    for (x = 0; x < l/2; ++x) {
        int y = buf[x]; buf[x] = buf[l-1-x]; buf[l-1-x] = y;
    }
    fputs(buf, fp);
    return 0;
}


char *cigar_from_bam(const bam1_t *b) {
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


/* Count matches and mismatches for an aligned read and also return
 * the corresponding qualities. returns NULL on error or pointer to
 * qual array for n_match and n_mismatch (sum is size). allocated
 * here. user must free
 */
int *count_matches(int *n_matches, int *n_mismatches,
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
int source_qual(const bam1_t *b, const char *ref)
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


static int mplp_func(void *data, bam1_t *b)
{
    extern int bam_realn(bam1_t *b, const char *ref);
    extern int bam_prob_realn_core(bam1_t *b, const char *ref, int);
    extern int bam_cap_mapQ(bam1_t *b, char *ref, int thres);
    mplp_aux_t *ma = (mplp_aux_t*)data;
    int ret, skip = 0;
    do {
        int has_ref;
        ret = ma->iter? bam_iter_read(ma->fp, ma->iter, b) : bam_read1(ma->fp, b);
        if (ret < 0) 
             break;
        if (b->core.tid < 0 || (b->core.flag&BAM_FUNMAP)) { /* exclude unmapped reads */
            skip = 1;
            continue;
        }
        if (ma->conf->bed) { /* test overlap */
            skip = !bed_overlap(ma->conf->bed, ma->h->target_name[b->core.tid], b->core.pos, bam_calend(&b->core, bam1_cigar(b)));
            if (skip) 
                 continue;
        }

        if (ma->conf->flag & MPLP_ILLUMINA13) {
            int i;
            uint8_t *qual = bam1_qual(b);
            for (i = 0; i < b->core.l_qseq; ++i)
                qual[i] = qual[i] > 31? qual[i] - 31 : 0;
        }
        has_ref = (ma->ref && ma->ref_id == b->core.tid)? 1 : 0;
        skip = 0;
        if (has_ref && (ma->conf->flag&MPLP_REALN)) 
             bam_prob_realn_core(b, ma->ref, (ma->conf->flag & MPLP_EXT_BAQ)? 3 : 1);
        if (has_ref && ma->conf->capQ_thres > 10) {
            int q = bam_cap_mapQ(b, ma->ref, ma->conf->capQ_thres);
            if (q < 0) {
                 skip = 1;
            } else if (b->core.qual > q) {
                 b->core.qual = q;
            }
        } else if (b->core.qual > ma->conf->max_mq) {
             b->core.qual = ma->conf->max_mq;
        } else if (b->core.qual < ma->conf->min_mq) {
             skip = 1; 
        }
        else if ((ma->conf->flag&MPLP_NO_ORPHAN) && (b->core.flag&1) && !(b->core.flag&2)) {
             skip = 1;
        }
    } while (skip);

    /* compute source qual if requested and have ref */
    if (ma->ref && ma->ref_id == b->core.tid && ma->conf->flag & MPLP_USE_SQ) {
         int sq = source_qual(b, ma->ref);
         LOG_FIXME("%s\n", "Got sq %d. What now?", sq);
    }

    return ret;
}



void process_plp(const bam_pileup1_t *plp, const int n_plp, 
                 const mplp_conf_t *conf, const char *ref, const int pos, 
                 const int ref_len, const char *target_name)
{
     int i;
     char ref_base;
     plp_col_t plp_col;

     /* "base counts" minus error-probs before base-level filtering
      * for each base. temporary data-structure for cheaply determining
      * consensus which is saved in plp_col */
     double base_counts[NUM_NT4] = { 0 }; 

     /* computation of depth (after read-level *and* base-level filtering)
      * samtools-0.1.18/bam2depth.c: 
      *   if (p->is_del || p->is_refskip) ++m; 
      *   else if (bam1_qual(p->b)[p->qpos] < bq) ++m
      * n_plp[i] - m
      */
     int num_skips = 0;

     ref_base = (ref && pos < ref_len)? ref[pos] : 'N';
     
     plp_col_init(& plp_col);
     plp_col.target = strdup(target_name);
     plp_col.pos = pos;
     plp_col.ref_base = toupper(ref_base);
     plp_col.coverage = n_plp;  /* this is coverage as in the original mpileup, 
                                   i.e. after read-level filtering */
     LOG_DEBUG("Processing %s:%d\n", plp_col.target, plp_col.pos+1);
     
     for (i = 0; i < n_plp; ++i) {
          /* Used parts of pileup_seq() here */
          const bam_pileup1_t *p = plp + i;
          int nt, nt4;
          int mq, bq; /* phred scores */
          int base_skip = 0; /* boolean */

          /* GATKs BI & BD: "are per-base quantities which estimate
           * the probability that the next base in the read was
           * mis-incorporated or mis-deleted (due to slippage, for
           * example)". See
           * http://www.broadinstitute.org/gatk/guide/article?id=44
           * and
           * http://gatkforums.broadinstitute.org/discussion/1619/baserecalibratorprintreads-bd-and-bi-flags
           */
          uint8_t *bi = bam_aux_get(p->b, "BI"); /* GATK indels */
          uint8_t *bd = bam_aux_get(p->b, "BD"); /* GATK deletions */

#if 0
          LOG_FIXME("At %s:%d %c: p->is_del=%d p->is_refskip=%d p->indel=%d p->is_head=%d p->is_tail=%d\n", 
                    plp_col.target, 
                    plp_col.pos+1,
                    bam_nt16_rev_table[bam1_seqi(bam1_seq(p->b), p->qpos)],
                    p->is_del, p->is_refskip, p->indel, p->is_head, p->is_tail);
#endif

          if (! p->is_del) {
               if (p->is_head) {
                    plp_col.num_heads += 1;
               }
               if (p->is_tail) {
                    plp_col.num_tails += 1;
               }

               /* nt for printing */
               nt = bam_nt16_rev_table[bam1_seqi(bam1_seq(p->b), p->qpos)];
               nt = bam1_strand(p->b)? tolower(nt) : toupper(nt);

               /* nt4 for indexing */
               nt4 = bam_nt16_nt4_table[bam1_seqi(bam1_seq(p->b), p->qpos)];                           

               bq = bam1_qual(p->b)[p->qpos];

               /* minimal base-call quality filtering
                */
               if (bq < conf->min_bq) {
                    base_skip = 1;
                    goto check_indel; /* FIXME: argh! */
               }

               /* the following will correct base-pairs down if they
                * exceed the valid sanger/phred limits. is it wise to
                * do this automatically? doesn't this indicate a
                * problem with the input ? */
               if (bq > 93) {
                    bq = 93; /* Sanger/Phred max */
               }

               base_counts[nt4] += (1.0 - PHREDQUAL_TO_PROB(bq));

               /* no need for check if mq is within user defined
                * limits. check was done in mplp_func */
               mq = p->b->core.qual;
               /* samtools check to detect Sanger max value: problem
                * is that an MQ Phred of 255 means NA according to the
                * samtools spec (needed below). This however is not
                * detectable if the following original samtools line
                * gets executed, which is why we remove it:
                * if (mq > 126) mq = 126;
                */

               PLP_COL_ADD_QUAL(& plp_col.base_quals[nt4], bq);
               PLP_COL_ADD_QUAL(& plp_col.map_quals[nt4], mq);
               if (bam1_strand(p->b)) {
                    plp_col.rv_counts[nt4] += 1;
               } else {
                    plp_col.fw_counts[nt4] += 1;
               }
                
               if (bi) {
                    int j;
                    char *t = (char*)(bi+1); /* 1 is type */
#if 0
                    printf("At %d qpos %d: BI=%s", plp_col.pos+1, p->qpos+1, t);
                    for (j = 0; j < p->indel; ++j) {
                         printf(" %c:%d-%d-%d", t[p->qpos+j], t[p->qpos+j-1]-33, t[p->qpos+j]-33, t[p->qpos+j+1]-33);
                    }
                    printf("\n");
#endif
                    /* only adding 1 value representing whole insert */
                    PLP_COL_ADD_QUAL(& plp_col.ins_quals, t[p->qpos]);
               }

               if (bd) {
                    int j;
                    char *t = (char*)(bd+1);  /* 1 is type */
#if 0
                    printf("At %d qpos %d: BD=%sx", plp_col.pos+1, p->qpos+1, t);
                    for (j = 0; j < p->indel; ++j) {
                         printf(" %c:%d-%d-%d", t[p->qpos+j], t[p->qpos+j-1]-33, t[p->qpos+j]-33, t[p->qpos+j+1]-33);
                    }
                    printf("\n");
#endif
                    PLP_COL_ADD_QUAL(& plp_col.del_quals, t[p->qpos]);
                    /* only adding 1 value representing whole insert */
               }
                            
          } /* ! p->is_del */
#if 0
          /* FIXME so what happens if is_del?
           * Observation: if indel<0 == del then they are followed by
           * is_del's. shouldn't that rather happen for ins? 
           */
          if (p->is_del || p->indel) {
               fflush(stdout); fflush(stderr);
               LOG_FIXME("At %s:%d p->is_del=%d p->is_refskip=%d p->indel=%d\n", 
                         plp_col.target, plp_col.pos+1, p->is_del, p->is_refskip, p->indel);
               fflush(stdout); fflush(stderr);

          }
#endif     
          
     check_indel:
          
          /* for post read- and base-level coverage */
          if (p->is_del || p->is_refskip || 1 == base_skip) {
               num_skips += 1;
          }
          
          /* A pattern \+[0-9]+[ACGTNacgtn]+' indicates there is an
           * insertion between this reference position and the next
           * reference position. The length of the insertion is given
           * by the integer in the pattern, followed by the inserted
           * sequence. Similarly, a pattern -[0-9]+[ACGTNacgtn]+
           * represents a deletion from the reference. The deleted
           * bases will be presented as ‘*’ in the following lines.
           */

          if (p->indel != 0) { /* seems to rule out is_del */
               int j;
               /* insertion (+)
                */
               if (p->indel > 0) {
                    plp_col.num_ins += 1;
                    plp_col.sum_ins += p->indel;
#ifdef PRINT_INDEL
                    putchar('+'); printw(p->indel, stdout);
                    putchar(' ');
                    for (j = 1; j <= p->indel; ++j) {
                         int c = bam_nt16_rev_table[bam1_seqi(bam1_seq(p->b), p->qpos + j)];
                         putchar(bam1_strand(p->b)? tolower(c) : toupper(c));
                    }
                    printf("\n");
#endif 
               /* deletion (-)
                */
               } else if (p->indel < 0) {
                    plp_col.num_dels += 1;
                    plp_col.sum_dels -= p->indel;
#ifdef PRINT_INDEL                   
                    printw(p->indel, stdout);
                    putchar(' ');
                    for (j = 1; j <= -p->indel; ++j) {
                         int c = (ref && (int)pos+j < ref_len)? ref[pos+j] : 'N';
                         putchar(bam1_strand(p->b)? tolower(c) : toupper(c));
                    }
                    printf("\n");
#endif
               }
          } /* if (p->indel != 0) ... */
          
     }  /* end: for (i = 0; i < n_plp; ++i) { */

     plp_col.coverage -= num_skips;
     
     /* determine consensus from 'counts' */
     plp_col.cons_base = bam_nt4_rev_table[
          argmax_d(base_counts, NUM_NT4)];

     if (debug) {
          plp_col_debug_print(& plp_col, stdout);
     }
#if 0
     plp_col_mpileup_print(& plp_col, conf, stdout);
#endif

     for (i = 0; i < NUM_NT4; ++i) {
          assert(plp_col.fw_counts[i] + plp_col.rv_counts[i] == plp_col.base_quals[i].n);
          assert(plp_col.base_quals[i].n == plp_col.map_quals[i].n);
          /* FIXME assert(plp_col.map_quals[i].n == plp_col.source_quals[i].n); */
     }

     call_lowfreq_snps(& plp_col, conf);

     plp_col_free(& plp_col);
}
/* process_plp */



static int mpileup(const mplp_conf_t *conf, const int n, const char **fn)
{
    mplp_aux_t **data;
    int i, tid, pos, *n_plp, tid0 = -1, beg0 = 0, end0 = 1u<<29, ref_len, ref_tid = -1, max_depth;
    const bam_pileup1_t **plp;
    bam_mplp_t iter;
    bam_header_t *h = 0;
    char *ref;
    kstring_t buf;
  
    /* paranoid exit. n only allowed to be one in our case (not much
     * of an *m*pileup, I know...) */
    if (1 != n) {
         fprintf(stderr, "FATAL(%s:%s): need exactly one BAM files as input (got %d)\n",
                 __FILE__, __FUNCTION__, n);
         for (i=0; i<n; i++) {
              fprintf(stderr, "%s\n", fn[i]);
         }
         exit(1);
    }

    memset(&buf, 0, sizeof(kstring_t));
    data = calloc(n, sizeof(void*));
    plp = calloc(n, sizeof(void*));
    n_plp = calloc(n, sizeof(int*));

#if 0
    fprintf(stderr, "[%s] Note: the format differs from regular pileup (see http://samtools.sourceforge.net/pileup.shtml) in the following ways\n", __func__);
    fprintf(stderr, "[%s] - bases and qualities are merged into one field\n", __func__);
    fprintf(stderr, "[%s] - each base is immediately followed by its quality\n", __func__);
    fprintf(stderr, "[%s] - on request mapping and base call quality are merged (P_joined = P_mq + (1-P_mq)*P_bq\n", __func__);
    fprintf(stderr, "[%s] - indel events are removed from bases and qualities and summarized in an additional field\n", __func__);
    fprintf(stderr, "[%s] - reference matches are not replaced with , or .\n", __func__);
#endif

    /* read the header and initialize data */
    for (i = 0; i < n; ++i) {
        bam_header_t *h_tmp;
        if (0 != strcmp(fn[i], "-")) {
          if (! file_exists(fn[i])) {
            fprintf(stderr, "File '%s' does not exist. Exiting...\n", fn[i]);
            exit(1);
          }
        }
        data[i] = calloc(1, sizeof(mplp_aux_t));
        data[i]->fp = strcmp(fn[i], "-") == 0? bam_dopen(fileno(stdin), "r") : bam_open(fn[i], "r");
        data[i]->conf = conf;
        h_tmp = bam_header_read(data[i]->fp);
        data[i]->h = i? h : h_tmp; /* for i==0, "h" has not been set yet */

        if (conf->reg) {
            int beg, end;
            bam_index_t *idx;
            idx = bam_index_load(fn[i]);
            if (idx == 0) {
                fprintf(stderr, "[%s] fail to load index for %d-th input.\n", __func__, i+1);
                exit(1);
            }
            if (bam_parse_region(h_tmp, conf->reg, &tid, &beg, &end) < 0) {
                fprintf(stderr, "[%s] malformatted region or wrong seqname for %d-th input.\n", __func__, i+1);
                exit(1);
            }
            if (i == 0) tid0 = tid, beg0 = beg, end0 = end;
            data[i]->iter = bam_iter_query(idx, tid, beg, end);
            bam_index_destroy(idx);
        }
        if (i == 0) {
             h = h_tmp;
        } else {
            bam_header_destroy(h_tmp);
        }
    }
    /*LOG_DEBUG("%s\n", "BAM header initialized");*/

    if (tid0 >= 0 && conf->fai) { /* region is set */
        ref = faidx_fetch_seq(conf->fai, h->target_name[tid0], 0, 0x7fffffff, &ref_len);
        ref_tid = tid0;
        for (i = 0; i < n; ++i) data[i]->ref = ref, data[i]->ref_id = tid0;
    } else {
         ref_tid = -1;
         ref = 0;
    }
    iter = bam_mplp_init(n, mplp_func, (void**)data);
    max_depth = conf->max_depth;
    if (max_depth * 1 > 1<<20)
        fprintf(stderr, "(%s) Max depth is above 1M. Potential memory hog!\n", __func__);
    if (max_depth * 1 < 8000) {
        max_depth = 8000 / 1;
        fprintf(stderr, "<%s> Set max per-file depth to %d\n", __func__, max_depth);
    }
    bam_mplp_set_maxcnt(iter, max_depth);

    vcf_write_header(conf->out, PACKAGE_STRING, conf->fa);

    LOG_DEBUG("%s\n", "Starting pileup loop");
    while (bam_mplp_auto(iter, &tid, &pos, n_plp, plp) > 0) {
         int i=0; /* NOTE: mpileup originally iterated over n */

        if (conf->reg && (pos < beg0 || pos >= end0))
             continue; /* out of the region requested */
        if (conf->bed && tid >= 0 && !bed_overlap(conf->bed, h->target_name[tid], pos, pos+1)) 
             continue;
        if (tid != ref_tid) {
            free(ref); ref = 0;
            if (conf->fai) 
                 ref = faidx_fetch_seq(conf->fai, h->target_name[tid], 0, 0x7fffffff, &ref_len);
            for (i = 0; i < n; ++i) 
                 data[i]->ref = ref, data[i]->ref_id = tid;
            ref_tid = tid;
        }
        i=0; /* i is 1 for first pos which is a bug due to the removal
              * of one of the loops, so reset here */

        process_plp(plp[i], n_plp[i], conf, 
                    ref, pos, ref_len, h->target_name[tid]);
    } /* while bam_mplp_auto */

    free(buf.s);
    bam_mplp_destroy(iter);
    bam_header_destroy(h);
    for (i = 0; i < n; ++i) {
        bam_close(data[i]->fp);
        if (data[i]->iter) bam_iter_destroy(data[i]->iter);
        free(data[i]);
    }
    free(data); free(plp); free(ref); free(n_plp);
    return 0;
}
/* mpileup */


void dump_mplp_conf(const mplp_conf_t *c, FILE *stream) {
     fprintf(stream, "mplp options\n");
     fprintf(stream, "  max_mq       = %d\n", c->max_mq);
     fprintf(stream, "  min_mq       = %d\n", c->min_mq);
     fprintf(stream, "  flag         = %d\n", c->flag);

     fprintf(stream, "  flag & MPLP_NO_ORPHAN  = %d\n", c->flag&MPLP_NO_ORPHAN?1:0);
     fprintf(stream, "  flag & MPLP_REALN      = %d\n", c->flag&MPLP_REALN?1:0);
     fprintf(stream, "  flag & MPLP_EXT_BAQ    = %d\n", c->flag&MPLP_EXT_BAQ?1:0);
     fprintf(stream, "  flag & MPLP_ILLUMINA13 = %d\n", c->flag&MPLP_ILLUMINA13?1:0);
     fprintf(stream, "  flag & MPLP_USE_MQ     = %d\n", c->flag&MPLP_USE_MQ?1:0);
     fprintf(stream, "  flag & MPLP_USE_SQ     = %d\n", c->flag&MPLP_USE_SQ?1:0);
     
     fprintf(stream, "  capQ_thres   = %d\n", c->capQ_thres);
     fprintf(stream, "  max_depth    = %d\n", c->max_depth);
     fprintf(stream, "  min_bq    = %d\n", c->min_bq);
     fprintf(stream, "  min_altbq = %d\n", c->min_altbq);
     fprintf(stream, "  def_altbq = %d\n", c->def_altbq);
     fprintf(stream, "  bonf         = %lu  (might get recalculated later)\n", c->bonf);
     fprintf(stream, "  sig          = %f\n", c->sig);
     fprintf(stream, "  reg          = %s\n", c->reg);
     fprintf(stream, "  fa           = %p\n", c->fa);
     /*fprintf(stream, "  fai          = %p\n", c->fai);*/
     fprintf(stream, "  bed          = %p\n", c->bed);
     fprintf(stream, "  cmdline      = %s\n", c->cmdline);
     fprintf(stream, "  out          = %p\n", c->out);
}


static void usage(const mplp_conf_t *mplp_conf)
{
     fprintf(stderr, "Usage: %s [options] in.bam\n\n", MYNAME);
     fprintf(stderr, "Options:\n");
     /* generic */
     fprintf(stderr, "          --verbose           be verbose\n");
     fprintf(stderr, "          --debug             enable debugging\n");
     /* regions */
     fprintf(stderr, "       -r|--region STR        region in which pileup should be generated [null]\n");
     fprintf(stderr, "       -l|--bed FILE          list of positions (chr pos) or regions (BED) [null]\n");
     /*  */
     fprintf(stderr, "       -d|--maxdepth INT      max per-BAM depth to avoid excessive memory usage [%d]\n", mplp_conf->max_depth);
     fprintf(stderr, "       -f|--reffa FILE        faidx indexed reference sequence file [null]\n");
     /* */
     fprintf(stderr, "       -o|--out FILE          vcf output file [- = stdout]\n");
     /* base call quality and baq */
     fprintf(stderr, "       -q|--min-bq INT        skip any base with baseQ smaller than INT [%d]\n", mplp_conf->min_bq);
     fprintf(stderr, "       -Q|--min-altbq INT     skip nonref-bases with baseQ smaller than INT [%d]. Not active if ref is N\n", mplp_conf->min_altbq);
     fprintf(stderr, "       -a|--def-altbq INT     nonref base qualities will be replace with this value [%d]\n", mplp_conf->def_altbq);
     fprintf(stderr, "       -B|--no-baq            disable BAQ computation\n");
     /* fprintf(stderr, "       -E           extended BAQ for higher sensitivity but lower specificity\n"); */
     /* mapping quality */
     fprintf(stderr, "       -m|--min_mq INT        skip alignments with mapQ smaller than INT [%d]\n", mplp_conf->min_mq);
     fprintf(stderr, "       -M|--max_mq INT        cap mapping quality at INT [%d]\n", mplp_conf->max_mq);
     fprintf(stderr, "       -J|--no-mq             don't merge mapQ into baseQ: P_e = P_mq + (1-P_mq) P_bq\n");
     fprintf(stderr, "       -S|--no-sq             don't merge sourceQ into baseQ\n");
     /* stats */
     fprintf(stderr, "       -s|--sig               P-value cutoff / significance level [%f]\n", mplp_conf->sig);
     fprintf(stderr, "       -b|--bonf              Bonferroni factor. INT or 'auto' (default; non-zero-cov-pos * 3)\n");
     fprintf(stderr, "                              'auto' needs to pre-parse the BAM file once, i.e. this won't work with input from stdin (or named pipes).\n");
     fprintf(stderr, "                              Higher numbers speed up computation on high-coverage data considerably.\n");
     /* misc */
     fprintf(stderr, "       -6|--illumina-1.3      assume the quality is Illumina-1.3-1.7/ASCII+64 encoded\n");
     fprintf(stderr, "       -A|--use-orphan        count anomalous read pairs\n");
}



int main_call(int argc, char *argv[])
{
     /* based on bam_mpileup() */
     int c, i;
     static int use_orphan = 0;
     int bonf_auto = 1;
     char *bam_file;
     char *bed_file = NULL;
     mplp_conf_t mplp;
     
     LOG_FIXME("%s\n", "- Proper source qual use missing");
     LOG_FIXME("%s\n", "- Indel handling missing");
     LOG_FIXME("%s\n", "- Implement routine test against old SNV caller");

    memset(&mplp, 0, sizeof(mplp_conf_t));
    /* default pileup options */
    mplp.max_mq = 255; /* 60 */
    mplp.min_bq = 3; /* 13 */
    mplp.min_altbq = 20; /* new */
    mplp.def_altbq = mplp.min_altbq;
    mplp.capQ_thres = 0;
    mplp.max_depth = 1000000; /* 250 */
    mplp.flag = MPLP_NO_ORPHAN | MPLP_REALN | MPLP_EXT_BAQ | MPLP_USE_MQ;/* | MPLP_USE_SQ; FIXME */
    mplp.bonf = 1;
    mplp.sig = 0.05;
    mplp.out = stdout;

    /* FIXME getopt should be replaced with something more sensible like
     * argtable2 ot Gopt. Otherwise there's always the risk between
     * incosistent long opt, short opt and usage */
    while (1) {
         static struct option long_opts[] = {
              /* see usage sync */
              {"verbose", no_argument, &verbose, 1},
              {"debug", no_argument, &debug, 1},

              {"region", required_argument, NULL, 'r'},
              {"bed", required_argument, NULL, 'l'},
              
              {"maxdepth", required_argument, NULL, 'd'},
              {"reffa", required_argument, NULL, 'f'},
               
              {"out", required_argument, NULL, 'o'},

              {"min-bq", required_argument, NULL, 'q'},
              {"min-altbq", required_argument, NULL, 'Q'},
              {"def-altbq", required_argument, NULL, 'a'},
              {"no-baq", no_argument, NULL, 'B'},
              /*{"ext-baq", required_argument, NULL, 'E'},*/
                   
              {"min-mq", required_argument, NULL, 'm'},
              {"max-mq", required_argument, NULL, 'M'},
              {"no-mq", no_argument, NULL, 'J'},
              {"no-sq", no_argument, NULL, 'S'},

              {"bonf", required_argument, NULL, 'b'},
              {"sig", required_argument, NULL, 's'},
                   
              {"illumina-1.3", no_argument, NULL, 'I'},
              {"use-orphan", no_argument, &use_orphan, 1},

              {"help", no_argument, NULL, 'h'},

              {0, 0, 0, 0} /* sentinel */
         };

         /* see usage sync */
         static const char *long_opts_str = "r:l:d:f:o:q:Q:a:Bm:M:JSb:s:IA:h";

         /* getopt_long stores the option index here. */
         int long_opts_index = 0;
         c = getopt_long(argc, argv, long_opts_str, long_opts, & long_opts_index);
         if (c == -1) {
              break;
         }

         switch (c) {
         /* see usage sync */
         case 'r': mplp.reg = strdup(optarg); break;
         case 'l': 
              mplp.bed = bed_read(optarg); 
              bed_file = strdup(optarg);
              break;
              
         case 'd': mplp.max_depth = atoi(optarg); break;
         case 'f':
              mplp.fa = strdup(optarg);
              mplp.fai = fai_load(optarg);
              if (mplp.fai == 0) 
                   return 1;
              break;
         case 'o':
              if (0 != strcmp(optarg, "-")) {
                   if (file_exists(optarg)) {
                        LOG_FATAL("Cowardly refusing to overwrite file '%s'. Exiting...\n", optarg);
                        return 1;
                   }
                   if (NULL == (mplp.out = fopen(optarg, "w"))) {
                        LOG_FATAL("Couldn't open file '%s'. Exiting...\n", optarg);
                        return 1;
                   }
              } /* else: already set to stdout */
         case 'q': mplp.min_bq = atoi(optarg); break;
         case 'Q': mplp.min_altbq = atoi(optarg); break;
         case 'a': mplp.def_altbq = atoi(optarg); break;
         case 'B': mplp.flag &= ~MPLP_REALN; break;
         /* case 'E': mplp.flag |= MPLP_EXT_BAQ; break; */
              
         case 'm': mplp.min_mq = atoi(optarg); break;
         case 'M': mplp.max_mq = atoi(optarg); break;
         case 'J': mplp.flag &= ~MPLP_USE_MQ; break;
         case 'S': mplp.flag &= ~MPLP_USE_SQ; break;
              
         case '6': mplp.flag |= MPLP_ILLUMINA13; break;

         case 'b': 
              if (0 == strncmp(optarg, "auto", 4)) {
                   bonf_auto = 1;

              } else {
                   bonf_auto = 0;
                   mplp.bonf = strtol(optarg, (char **)NULL, 10); /* atol */ 
                   if (0==mplp.bonf) {
                        LOG_FATAL("%s\n", "Couldn't parse Bonferroni factor\n"); 
                        exit(1);
                   }
              }
              break;
         case 's': 
              mplp.sig = strtof(optarg, (char **)NULL); /* atof */
              if (0==mplp.sig) {
                   LOG_FATAL("%s\n", "Couldn't parse sign-threshold\n"); 
                   exit(1);
              }
              break;
              
         case 'h': usage(& mplp); exit(0); /* WARN: not printing defaults if some args where parsed */
         case '?': LOG_FATAL("%s\n", "unrecognized arguments found. Exiting...\n"); exit(1);
#if 0
         case 0:  fprintf(stderr, "ERROR: long opt (%s) not mapping to short option. Exiting...\n", long_opts[long_opts_index].name); exit(1);
#endif
         default:
              break;
         }
    }
    if (use_orphan) {
         mplp.flag &= ~MPLP_NO_ORPHAN;
    }
    mplp.cmdline[0] = '\0';
    for (i=0; i<argc; i++) {
         strncat(mplp.cmdline, argv[i], sizeof(mplp.cmdline)-strlen(mplp.cmdline)-2);
         strcat(mplp.cmdline, " ");
    }

    if (argc == 1) {
        fprintf(stderr, "\n");
        usage(& mplp);
        return 1;
    }
    if (1 != argc - optind) {
        fprintf(stderr, "Need exactly one BAM file as last argument\n");
        return(EXIT_FAILURE);
    }
    bam_file = (argv + optind)[0];


    if (bonf_auto) {
         double cov_mean;
         long int num_non0cov_pos;
         LOG_DEBUG("Automatically determining Bonferroni factor for bam=%s reg=%s bed=%s\n",
                   bam_file, mplp.reg, bed_file); 
         if (depth_stats(&cov_mean, &num_non0cov_pos, bam_file, mplp.reg, bed_file,
                         &mplp.min_bq, &mplp.min_mq)) {
              LOG_FATAL("%s\n", "Couldn't determine Bonferroni factor automatically\n"); 
              exit(1);
         }
         mplp.bonf = num_non0cov_pos*3;
         LOG_VERBOSE("Automatically determined Bonferroni factor = %lu\n", mplp.bonf);
    }


    if (debug) {
         dump_mplp_conf(& mplp, stderr);
    }

    /* FIXME: implement logic_check_opts() */
    assert(mplp.min_mq <= mplp.max_mq);
    assert(mplp.min_bq <= mplp.min_altbq);
    assert(! (mplp.bed && mplp.reg));
   
    mpileup(&mplp, 1, (const char **) argv + optind);

    if (mplp.out == stdout) {
         fclose(mplp.out);
    }
    free(mplp.reg); 
    free(mplp.fa);
    if (mplp.fai) {
         fai_destroy(mplp.fai);
    }
    free(bed_file);
    if (mplp.bed) {
         bed_destroy(mplp.bed);
    }
    LOG_VERBOSE("%s\n", "Successful exit.");
    return 0;
}
/* main_call */


#ifdef MAIN_TEST

int main()
{
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

