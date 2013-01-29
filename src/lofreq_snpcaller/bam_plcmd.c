/* -*- c-file-style: "k&r" -*-
 * http://emacswiki.org/emacs/IndentingC
 * http://www.emacswiki.org/emacs/LocalVariables
 * http://en.wikipedia.org/wiki/Indent_style 
 *
 * This is based on samtools bam_plcmd.c
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
#include <getopt.h>

#include "sam.h"
#include "faidx.h"
#include "kstring.h"
#include "snpcaller.h"
#include "bam2depth.h"
#include "utils.h">
/* from bedidx.c */
void *bed_read(const char *fn);
void bed_destroy(void *_h);
int bed_overlap(const void *_h, const char *chr, int beg, int end);


/* mpileup configuration flags 
 */
/*#define MPLP_GLF   0x10*/
#define MPLP_NO_COMP 0x20
#define MPLP_NO_ORPHAN 0x40
#define MPLP_REALN   0x80
#define MPLP_FMT_DP 0x100
#define MPLP_FMT_SP 0x200
#define MPLP_NO_INDEL 0x400
#define MPLP_EXT_BAQ 0x800
#define MPLP_ILLUMINA13 0x1000
/*#define MPLP_IGNORE_RG 0x2000
#define MPLP_PRINT_POS 0x4000*/
#define MPLP_PRINT_MAPQ 0x8000
#define MPLP_JOIN_BQ_AND_MQ 0x10000


#define MYNAME "lofreq_mpileup"
#define PHREDQUAL_TO_PROB(phred) (pow(10.0, -1.0*(phred)/10.0))
#define PROB_TO_PHREDQUAL(prob) ((int)(-10.0 * log10(prob)))


const char *bam_nt4_rev_table = "ACGTN"; /* similar to bam_nt16_rev_table */
#define NUM_NT4 5 /* strlen(bam_nt4_rev_table); */



/* mpileup configuration structure 
 */
typedef struct {
     int max_mq, min_mq;
     int flag;
     int capQ_thres;
     int max_depth;
     int min_baseQ, min_altbaseQ;
     int def_altbaseQ;
     unsigned long int bonf;
     double sig;
     char *reg;
     faidx_t *fai;
     void *bed;
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



/* Logging macros
 *
 * Taken from squicl-0.2.8
 * You must use at least one fmt+string and append trailing "\n", e.g.
 * "%s\n", "string", instead of "string\n"
 *
 */
int debug = 0;
int verbose = 0;
/* Taken from the Linux kernel source and slightly modified.
 * bool_flag: print or don't
 */
int
printk(FILE *stream, int bool_flag, const char *fmt, ...)
{                
    va_list args;
    static char printk_buf[8192];
    int printed_len=0;
 
    if (bool_flag) {
        /* Emit the output into the temporary buffer */
        va_start(args, fmt);
        printed_len = vsnprintf(printk_buf, sizeof(printk_buf), fmt, args);
        va_end(args);

        fprintf(stream, "%s", printk_buf);
        fflush(stream);        
    }
    return printed_len;
}
/* print only if debug is true*/
#define LOG_DEBUG(fmt, args...)     printk(stderr, debug, "DEBUG(%s|%s): " fmt, __FILE__, __FUNCTION__, ## args)
/* print only if verbose is true*/
#define LOG_VERBOSE(fmt, args...)   printk(stderr, verbose || debug, fmt, ## args)
/* always warn to stderr */
#define LOG_WARN(fmt, args...)      printk(stderr, 1, "WARNING(%s|%s): " fmt, __FILE__, __FUNCTION__, ## args)
/* always print errors to stderr*/
#define LOG_ERROR(fmt, args...)     printk(stderr, 1, "ERROR(%s|%s): " fmt, __FILE__, __FUNCTION__, ## args)
/* always print critical errors to stderr*/
#define LOG_CRITICAL(fmt, args...)  printk(stderr, 1, "CRITICAL(%s|%s): " fmt, __FILE__, __FUNCTION__, ## args)
#define LOG_FATAL(fmt, args...)     printk(stderr, 1, "FATAL(%s|%s): " fmt, __FILE__, __FUNCTION__, ## args)
/* always print fixme's */
#define LOG_FIXME(fmt, args...)  printk(stderr, 1, "FIXME(%s|%s:%d): " fmt, __FILE__, __FUNCTION__, __LINE__, ## args)


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

     /* FIXME only temporary before they move into they own structure */
     int num_ins, sum_ins;
     int num_dels, sum_dels;
} plp_col_t;


#define PLP_COL_ADD_QUAL(p, q)   int_varray_add_value((p), (q));


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
         int_varray_init(& p->source_quals[i], 0);
         p->fw_counts[i] = 0;
         p->rv_counts[i] = 0;
    }

    p->num_heads = p->num_tails = 0;

    p->num_ins = p->sum_ins = 0;
    p->num_dels = p->sum_dels = 0;

}

void plp_col_free(plp_col_t *p) {
    int i;

    free(p->target);
    for (i=0; i<NUM_NT4; i++) {
         int_varray_free(& p->base_quals[i]);
         int_varray_free(& p->map_quals[i]);
         int_varray_free(& p->source_quals[i]);
    }
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
void plp_col_mpileup_print(const plp_col_t *p, mplp_conf_t *conf, FILE *stream)
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


/* "Merge" MQ and BQ if requested and if MAQP not 255 (not available):
 *  P_jq = P_mq * + (1-P_mq) P_bq.
 */
int merge_baseq_and_mapq(int bq, int mq)
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
 * Assuming conf->min_baseQ and read-level filtering was already done
 * upstream. altbase mangling happens here however.
 * 
 */
void call_lowfreq_snps(const plp_col_t *p, mplp_conf_t *conf)
{
     int *quals; /* qualities passed down to snpcaller */
     int quals_len; /* #elements in quals */
     int i, j;

     /* 4 bases ignoring N, -1 reference/consensus base makes 3 */
     double pvalues[3]; /* pvalues reported back from snpcaller */
     int alt_counts[3]; /* counts for alt bases handed down to snpcaller */
     int alt_raw_counts[3];
     int alt_bases[3];/* actual alt bases */
     int alt_idx;
     int got_alt_bases = 0;

     /* don't call if no coverage or if we don't know what to call
      * against */
     if (p->coverage == 0 || p->cons_base == 'N') {          
          return;
     }

     /* check for consensus snps, i.e. those where the consensus
      * determined here is different from the reference coming from a
      * fasta file */
     if (p->ref_base != 'N' && p->ref_base != p->cons_base) {
          LOG_FIXME("cons var snp: %s %d %c>%c\n", p->target, p->pos+1, p->ref_base, p->cons_base);
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
                    if (bq < conf->min_altbaseQ) {
                         continue; /* WARNING base counts now invalid. We used them for freq reporting anyway, otherwise heterozygous calls look odd */
                    }
                    bq = conf->def_altbaseQ;
                    alt_counts[alt_idx] += 1;
               }

               if ((conf->flag & MPLP_JOIN_BQ_AND_MQ)) {
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
               LOG_FIXME("low freq snp: %s %d %c>%c pv:%f;raw:%d/%d;filt:%d/%d\n",
                         p->target, p->pos+1, p->cons_base, alt_base,
                         pvalue, 
                         alt_raw_count, p->coverage,
                         alt_count, quals_len);
          }
     }
     free(quals);
}



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
    return ret;
}



void process_plp(const bam_pileup1_t *plp, const int n_plp, 
                 mplp_conf_t *conf, const char *ref, const int pos, 
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
      *   else if (bam1_qual(p->b)[p->qpos] < baseQ) ++m
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

     for (i = 0; i < n_plp; ++i) {
          /* used parts of pileup_seq() here */
          const bam_pileup1_t *p = plp + i;
          int nt, nt4;
          int mq, bq; /* phred scores */
          int base_skip = 0; /* boolean */
             
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
               if (bq < conf->min_baseQ) {
                    base_skip = 1;
                    goto check_indel; /* FIXME: argh! */
               }

               /* the following will correct base-pairs down if they
                * exceed the valid sanger/phred limits. is it wise to
                * do this automatically? doesn't this indicate a
                * problem with the input ? */
               if (bq > 93) 
                    bq = 93; /* Sanger/Phred max */

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
                                            
          } /* ! p->is_del */
          
          
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
          if (p->indel != 0) {
               if (p->indel > 0) {
                    /* + */
                    plp_col.num_ins += 1;
                    plp_col.sum_ins += p->indel;
                    
#ifdef INDEL_TEST
                    {int j;
                         uint8_t *s = bam_aux_get(p->b, "BI");
                         printf("\n");
                         if (s) {
                              char *t = (char*)(s+1);
                              printf("BI:%s ", t);
                              for (j = 0; j < p->indel; ++j) {
                                   printf("%c", t[p->qpos+j]);
                              }
                              printf("\n");
                         }
                         
                         putchar('+'); printw(p->indel, stdout);
                         putchar(' ');
                         for (j = 1; j <= p->indel; ++j) {
                              int c = bam_nt16_rev_table[bam1_seqi(bam1_seq(p->b), p->qpos + j)];
                              putchar(bam1_strand(p->b)? tolower(c) : toupper(c));
                         }
                         printf("\n");
                    }
#endif
                    
               } else if (p->indel < 0) {
                    /* - */
                    plp_col.num_dels += 1;
                    plp_col.sum_dels -= p->indel;
                    
#ifdef INDEL_TEST
                    {int j;
                         uint8_t *s = bam_aux_get(p->b, "BD");
                         printf("\n");
                         if (s) {
                              char *t = (char*)(s+1);
                              printf("BD:%s ", t);
                              for (j = 0; j < -p->indel; ++j) {
                                   printf("%c", t[p->qpos+j]);
                              }
                              printf("\n");
                         }
                         
                         printw(p->indel, stdout);
                         putchar(' ');
                         for (j = 1; j <= -p->indel; ++j) {
                              int c = (ref && (int)pos+j < ref_len)? ref[pos+j] : 'N';
                              putchar(bam1_strand(p->b)? tolower(c) : toupper(c));
                         }
                         printf("\n");
                    }
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



static int mpileup(mplp_conf_t *conf, int n, char **fn)
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
        if (i == 0) h = h_tmp;
        else {
             /* FIXME: to check consistency */
            bam_header_destroy(h_tmp);
        }
    }

    if (tid0 >= 0 && conf->fai) { /* region is set */
        ref = faidx_fetch_seq(conf->fai, h->target_name[tid0], 0, 0x7fffffff, &ref_len);
        ref_tid = tid0;
        for (i = 0; i < n; ++i) data[i]->ref = ref, data[i]->ref_id = tid0;
    } else ref_tid = -1, ref = 0;
    iter = bam_mplp_init(n, mplp_func, (void**)data);
    max_depth = conf->max_depth;
    if (max_depth * 1 > 1<<20)
        fprintf(stderr, "(%s) Max depth is above 1M. Potential memory hog!\n", __func__);
    if (max_depth * 1 < 8000) {
        max_depth = 8000 / 1;
        fprintf(stderr, "<%s> Set max per-file depth to %d\n", __func__, max_depth);
    }
    bam_mplp_set_maxcnt(iter, max_depth);


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



void dump_mplp_conf(const mplp_conf_t *c, FILE *stream) {
     fprintf(stream, "mplp options\n");
     fprintf(stream, " max_mq       = %d\n", c->max_mq);
     fprintf(stream, " min_mq       = %d\n", c->min_mq);
     fprintf(stream, " flag         = %d\n", c->flag);
     fprintf(stream, " capQ_thres   = %d\n", c->capQ_thres);
     fprintf(stream, " max_depth    = %d\n", c->max_depth);
     fprintf(stream, " min_baseQ    = %d\n", c->min_baseQ);
     fprintf(stream, " min_altbaseQ = %d\n", c->min_altbaseQ);
     fprintf(stream, " def_altbaseQ = %d\n", c->def_altbaseQ);
     fprintf(stream, " bonf         = %lu\n", c->bonf);
     fprintf(stream, " sig          = %f\n", c->sig);
     fprintf(stream, " reg          = %s\n", c->reg);
     fprintf(stream, " fai          = %p\n", c->fai);
     fprintf(stream, " bed          = %p\n", c->bed);
}


void usage(const mplp_conf_t *mplp_conf) {
     fprintf(stderr, "Usage: %s [mpileup] [options] in.bam\n\n", MYNAME);
     fprintf(stderr, "Options:\n");
     /* generic */
     fprintf(stderr, "          --verbose           be verbose\n");
     fprintf(stderr, "          --debug             enable debugging\n");
     /* regions */
     fprintf(stderr, "       -r|--region STR        region in which pileup is generated [null]\n");
     fprintf(stderr, "       -l|--bed FILE          list of positions (chr pos) or regions (BED) [null]\n");
     /*  */
     fprintf(stderr, "       -d|--maxdepth INT      max per-BAM depth to avoid excessive memory usage [%d]\n", mplp_conf->max_depth);
     fprintf(stderr, "       -f|--reffa FILE        faidx indexed reference sequence file [null]\n");
     /* base call quality and baq */
     fprintf(stderr, "       -q|--min_baseq INT     skip any base with baseQ smaller than INT [%d]\n", mplp_conf->min_baseQ);
     fprintf(stderr, "       -Q|--min_altbaseq INT  skip nonref-bases with baseQ smaller than INT [%d]. Not active if ref is N\n", mplp_conf->min_altbaseQ);
     fprintf(stderr, "       -a|--def_altbaseq INT  nonref base qualities will be replace with this value [%d]\n", mplp_conf->def_altbaseQ);
     fprintf(stderr, "       -B|--no-baq            disable BAQ computation\n");
     /* fprintf(stderr, "       -E           extended BAQ for higher sensitivity but lower specificity\n"); */
     /* mapping quality */
     fprintf(stderr, "       -m|--min_mq INT        skip alignments with mapQ smaller than INT [%d]\n", mplp_conf->min_mq);
     fprintf(stderr, "       -M|--max_mq INT        cap mapping quality at INT [%d]\n", mplp_conf->max_mq);
     fprintf(stderr, "       -j|--join-quals        join mapQ and baseQ per base: P_e = P_mq + (1-P_mq) P_bq\n");
     /* stats */
     fprintf(stderr, "       -s|--sig               P-value cutoff / significance level [%f]\n", mplp_conf->sig);
     fprintf(stderr, "       -b|--bonf              Bonferroni factor [%lu]. INT or 'auto' (non-zero-cov-pos * 3\n", mplp_conf->bonf);
     fprintf(stderr, "                              'auto' needs to pre-parse BAM, i.e. won't work with input from stdin.\n");
     fprintf(stderr, "                              Higher numbers speed up computation on high-coverage data considerably.\n");
     /* misc */
     fprintf(stderr, "       -6|--illumina-1.3      assume the quality is Illumina-1.3-1.7/ASCII+64 encoded\n");
     fprintf(stderr, "       -A|--use-orphan        count anomalous read pairs\n");

     fprintf(stderr, "\nDefault parameters here differ from original samtools mpileup\n");
     fprintf(stderr, "Furthermore, the format used here differs from regular pileup (see http://samtools.sourceforge.net/pileup.shtml):\n");
     fprintf(stderr, " - bases and qualities are merged into one field and each base is immediately followed by its quality\n");
     fprintf(stderr, " - on request mapping and base call quality are merged (P_joined = P_mq + (1-P_mq)*P_bq\n");
     fprintf(stderr, " - indel events are removed from bases and qualities and summarized in an additional field\n");
     fprintf(stderr, " - reference matches are not replaced with , or .\n\n");
}



int bam_mpileup(int argc, char *argv[])
{
    int c;
    static int use_orphan = 0;
    int bonf_auto = 0;
    char *bam_file;
    char *bed_file = NULL;
    mplp_conf_t mplp;

    memset(&mplp, 0, sizeof(mplp_conf_t));
    /* default pileup options */
    mplp.max_mq = 255; /* 60 */
    mplp.min_baseQ = 3; /* 13 */
    mplp.min_altbaseQ = 20; /* new */
    mplp.def_altbaseQ = mplp.min_altbaseQ;
    mplp.capQ_thres = 0;
    mplp.max_depth = 1000000; /* 250 */
    /* mplp.flag = MPLP_NO_ORPHAN | MPLP_REALN; */
    mplp.flag = MPLP_NO_ORPHAN | MPLP_REALN | MPLP_EXT_BAQ;
    mplp.bonf = 1;
    mplp.sig = 0.05;
    /* should differentiate between pileup and snp calling options */

    while (1) {
         static struct option long_opts[] = {
              /* These options set a flag. */
              {"verbose", no_argument, &verbose, 1},
              {"debug", no_argument, &debug, 1},

              {"region", required_argument, NULL, 'r'},
              {"bed", required_argument, NULL, 'l'},
              
              {"maxdepth", required_argument, NULL, 'd'},
              {"reffa", required_argument, NULL, 'f'},
               
              {"min_baseq", required_argument, NULL, 'q'},
              {"min_altbaseq", required_argument, NULL, 'Q'},
              {"def_altbaseq", required_argument, NULL, 0},
              {"no-baq", no_argument, NULL, 'B'},
              /*{"ext-baq", required_argument, NULL, 'E'},*/
                   
              {"min_mq", required_argument, NULL, 'm'},
              {"max_mq", required_argument, NULL, 'M'},
              {"join_quals", no_argument, NULL, 'j'},

              {"bonf", required_argument, NULL, 'b'},
              {"sig", required_argument, NULL, 's'},
                   
              {"illumina-1.3", no_argument, NULL, '6'},
              {"use-orphan", no_argument, &use_orphan, 1},

              {0, 0, 0, 0} /* sentinel */
         };
         /* WARN keep in sync with above */
         static const char *long_opts_str = "r:l:d:f:Q:q:Bm:M:jb:s:6A:h"; 

         /* getopt_long stores the option index here. */
         int long_opts_index = 0;
         c = getopt_long(argc, argv, long_opts_str, long_opts, & long_opts_index);
         if (c == -1) {
              break;
         }

         switch (c) {
         case 'r': mplp.reg = strdup(optarg); break;
         case 'l': 
              mplp.bed = bed_read(optarg); 
              bed_file = strdup(optarg);
              break;
              
         case 'd': mplp.max_depth = atoi(optarg); break;
         case 'f':
              mplp.fai = fai_load(optarg);
              if (mplp.fai == 0) 
                   return 1;
              break;
              
         case 'q': mplp.min_baseQ = atoi(optarg); break;
         case 'Q': mplp.min_altbaseQ = atoi(optarg); break;
         case 'a': mplp.def_altbaseQ = atoi(optarg); break;
         case 'B': mplp.flag &= ~MPLP_REALN; break;
         /* case 'E': mplp.flag |= MPLP_EXT_BAQ; break; */
              
         case 'm': mplp.min_mq = atoi(optarg); break;
         case 'M': mplp.max_mq = atoi(optarg); break;
         case 'j': mplp.flag |= MPLP_JOIN_BQ_AND_MQ; break;
              
         case '6': mplp.flag |= MPLP_ILLUMINA13; break;

         case 'b': 
              if (0 == strncmp(optarg, "auto", 4)) {
                   bonf_auto = 1;

              } else {
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
    LOG_FIXME("bam_file = %s\n", bam_file);


    if (bonf_auto) {
         double cov_mean;
         long int num_non0cov_pos;
         LOG_DEBUG("Automatically determining Bonferroni factor for bam=%s reg=%s bed=%s\n",
                   bam_file, mplp.reg, bed_file); 
         if (depth_stats(&cov_mean, &num_non0cov_pos, bam_file, mplp.reg, bed_file,
                         &mplp.min_baseQ, &mplp.min_mq)) {
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
    assert(mplp.min_baseQ <= mplp.min_altbaseQ);
    assert(! (mplp.bed && mplp.reg));
   


    mpileup(&mplp, 1, argv + optind);

    free(mplp.reg); 
    if (mplp.fai) {
         fai_destroy(mplp.fai);
    }
    free(bed_file);
    if (mplp.bed) {
         bed_destroy(mplp.bed);
    }

    return 0;
}



int main(int argc, char **argv)
{
     LOG_FIXME("%s\n", "- Source qual missing");
     LOG_FIXME("%s\n", "- Missing test against old SNV caller");
     LOG_FIXME("%s\n", "- Set defaults once things are working");

     if (argc>1 && strcmp(argv[1], "mpileup") == 0) {
          return bam_mpileup(argc-1, argv+1);
     } else {
          return bam_mpileup(argc, argv);
     }
}



#ifdef SANDBOX
void sandbox()
{
     char test_nucs[] = "ACGTNRYacgtnryZ\0";
     int i;
     
     for (i=0; i<strlen(test_nucs); i++) {
          printf("%d %c - %d - %d\n", i, test_nucs[i],
                 bam_nt16_table[(int)test_nucs[i]],
                 bam_nt16_nt4_table[bam_nt16_table[(int)test_nucs[i]]]);
     }
}
#endif

