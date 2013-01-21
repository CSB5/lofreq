/* -*- c-file-style: "k&r" -*-
 * http://emacswiki.org/emacs/IndentingC
 * http://www.emacswiki.org/emacs/LocalVariables
 * http://en.wikipedia.org/wiki/Indent_style 
 *
 * This is a modified version of samtools mpileup, based on
 * bam_plcmd.c.
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

#include "sam.h"
#include "faidx.h"
#include "kstring.h"
/* from bedidx.c */
void *bed_read(const char *fn);
void bed_destroy(void *_h);
int bed_overlap(const void *_h, const char *chr, int beg, int end);



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


typedef struct {
    int max_mq, min_mq, flag, capQ_thres, max_depth;
    int min_baseQ, min_altbaseQ; /* new */
    char *reg, *pl_list;
    faidx_t *fai;
    void *bed, *rghash;
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



/* ------------------------------ */

#define MYNAME "lofreq_mpileup"
#define PHREDQUAL_TO_PROB(phred) (pow(10.0, -1.0*(phred)/10.0))
#define PROB_TO_PHREDQUAL(prob) ((int)(-10.0 * log10(prob)))


const char *bam_nt4_rev_table = "ACGTN"; /* as bam_nt16_rev_table */
#define NUM_NT4 5 /* strlen(bam_nt4_rev_table); */

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


/* Logging macros
 *
 * Taken from squicl-0.2.8
 * You must use at least one fmt+string and append trailing "\n", e.g.
 * "%s\n", "string", instead of "string\n"
 *
 */
int debug = 0;
int verbose = 1;
/* print only if debug is true*/
#define LOG_DEBUG(fmt, args...)     printk(stderr, debug, "DEBUG(%s|%s): " fmt, __FILE__, __FUNCTION__, ## args)
/* print only if verbose is true*/
#define LOG_VERBOSE(fmt, args...)   printk(stderr, verbose, fmt, ## args)
/* always warn to stderr */
#define LOG_WARN(fmt, args...)      printk(stderr, 1, "WARNING(%s|%s): " fmt, __FILE__, __FUNCTION__, ## args)
/* always print errors to stderr*/
#define LOG_ERROR(fmt, args...)     printk(stderr, 1, "ERROR(%s|%s): " fmt, __FILE__, __FUNCTION__, ## args)
/* always print critical errors to stderr*/
#define LOG_CRITICAL(fmt, args...)     printk(stderr, 1, "CRITICAL(%s|%s): " fmt, __FILE__, __FUNCTION__, ## args)
/* always print fixme's */
#define LOG_FIXME(fmt, args...)  printk(stderr, 1, "FIXME(%s|%s:%d): " fmt, __FILE__, __FUNCTION__, __LINE__, ## args)


typedef struct {
     unsigned long int n; /* number of elements stored */
     int *data; /* actual array of data */

     size_t grow_by_size; /* if needed grow array by this value. will double previous size if <=1 */
     size_t alloced; /* actually allocated size for data */
} int_varray_t;



void int_varray_free(int_varray_t *a) 
{
    assert(NULL != a);

    free(a->data); /* save even if a->data==NULL */
    a->data = NULL;
    a->n = a->alloced = a->grow_by_size = 0;
}

void int_varray_init(int_varray_t *a, 
                     const size_t grow_by_size)
{
    assert(NULL != a);
    assert(0 == a->alloced); /* otherwise use clear() */

    a->n = 0;
    a->data = NULL;

    a->grow_by_size = grow_by_size;
    a->alloced = 0;
}

void int_varray_add_value(int_varray_t *a, const int value)
{
    assert(NULL != a);
    if (a->n * sizeof(int) == a->alloced) {
        size_t size_to_alloc;
        
        if (1 >=  a->grow_by_size) {
             assert(SIZE_MAX - a->alloced > a->alloced);
             size_to_alloc = 0==a->n ? sizeof(int) : a->alloced*2;
        } else {
             assert(SIZE_MAX - a->alloced > a->grow_by_size);
             size_to_alloc = a->alloced + a->grow_by_size;
        }
#ifdef DEBUG
        LOG_DEBUG("size=%lu alloced=%zd size_to_alloc=%zd sizeof(data)=%zd\n",
                a->n, a->alloced, size_to_alloc, sizeof(a->data));
#endif
        a->data = realloc(a->data, size_to_alloc);
        a->alloced = size_to_alloc;
    }

    a->data[a->n] = value;
    a->n++;
}


typedef struct {
     char *target; /* chromsome or sequence name */
     int pos; /* position */
     char ref_base; /* reference base */
     int coverage; /* coverage after read-level filtering (same as in samtools mpileup (n_plp)) */
     int_varray_t qs_fw[NUM_NT4]; /* qualities on fw strand for all bases */
     int_varray_t qs_rv[NUM_NT4]; /* qualities on rv strand for all bases */
     int num_heads; /* number of read starts at this pos */
     int num_tails; /* number of read ends at this pos */

     /* FIXME only temporary before they move into they own structure */
     int num_ins, sum_ins;
     int num_dels, sum_dels;

} plp_col_t;


void plp_col_init(plp_col_t *p) {
    int i;

    p->target =  NULL;
    p->pos = -INT_MAX;
    p->ref_base = '\0';
    p->coverage = -INT_MAX;
    for (i=0; i<NUM_NT4; i++) {
         int_varray_init(& p->qs_fw[i], 0);
         int_varray_init(& p->qs_rv[i], 0);
    }
    p->num_heads = p->num_tails = 0;

    p->num_ins = p->sum_ins = 0;
    p->num_dels = p->sum_dels = 0;

}

void plp_col_free(plp_col_t *p) {
    int i;

    free(p->target);
    for (i=0; i<NUM_NT4; i++) {
         int_varray_free(& p->qs_fw[i]);
         int_varray_free(& p->qs_rv[i]);
    }
}

void plp_col_debug_print(const plp_col_t *p, FILE *stream) {
     int i;
     
     fprintf(stream, "%s\t%d\t%c base-counts (fw/rv): ", 
             p->target, p->pos+1, p->ref_base);
     for (i=0; i<NUM_NT4; i++) {
          fprintf(stream, " %c:%lu/%lu",
                  bam_nt4_rev_table[i],
                  p->qs_fw[i].n,
                  p->qs_rv[i].n);
     }
     fprintf(stream, " heads:%d tails:%d\n", p->num_heads, p->num_tails);
     fprintf(stream, "\n");

}

void plp_col_mpileup_print(const plp_col_t *p, FILE *stream) {
     int i;
     
     fprintf(stream, "%s\t%d\t%c\t%d\t",
             p->target, p->pos+1, p->ref_base, p->coverage);

     for (i=0; i<NUM_NT4; i++) {
          int j, nt, q;

          nt = toupper(bam_nt4_rev_table[i]);
          for (j=0; j<p->qs_fw[i].n; j++) {
               q = p->qs_fw[i].data[j];
               fprintf(stream, "%c%c", nt, q+33);
          }

          nt = tolower(bam_nt4_rev_table[i]);
          for (j=0; j<p->qs_rv[i].n; j++) {
               q = p->qs_rv[i].data[j];
               fprintf(stream, "%c%c", nt, q+33);
          }
     }             

     fprintf(stream, "\t#heads=%d #tails=%d #ins=%d ins_len=%.1f #del=%d del_len=%.1f",
            p->num_heads, p->num_tails,
            p->num_ins, p->num_ins ? p->sum_ins/(float)p->num_ins : 0,
            p->num_dels, p->num_dels ? p->sum_dels/(float)p->num_dels : 0);

     fprintf(stream, "\n");
}

/* p should be a pointer to the quals of a nuc of a certrain strand, e.g.
 * & plp_col.qs_rv[nt4]. q is the quality.
*/
#define PLP_COL_ADD_BASEQUAL(p, q)   int_varray_add_value((p), (q));



/* FIXME also in lofreq_core:utils.c */
static inline int file_exists(char *fname) 
{
  /* from 
   * http://stackoverflow.com/questions/230062/whats-the-best-way-to-check-if-a-file-exists-in-c-cross-platform 
   */
  if (access(fname, F_OK) != -1) {
      return 1;
  } else {
      return 0;
  }
}

/* ------------------------------ */



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
     int j;
     char ref_base;
     plp_col_t plp_col; /* FIXME in the long-run meant to be the main data-structure */

     /* computation of depth (after read-level *and* base-level filtering)
      * samtools-0.1.18/bam2depth.c: 
      *   if (p->is_del || p->is_refskip) ++m; 
      *   else if (bam1_qual(p->b)[p->qpos] < baseQ) ++m
      * n_plp[i] - m
      */
     int f_depth = 0;
     int num_skips = 0;

     ref_base = (ref && pos < ref_len)? ref[pos] : 'N';
     
     plp_col_init(& plp_col);
     plp_col.target = strdup(target_name);
     plp_col.pos = pos;
     plp_col.ref_base = ref_base;
     plp_col.coverage = n_plp;  /* this is a in mpileup the coverage, but after read-level filtering. */

     for (j = 0; j < n_plp; ++j) {
          /* used parts of pileup_seq() here */
          const bam_pileup1_t *p = plp + j;
          int nt, nt4;
          int mq, bq, jq, final_q; /* phred scores */
          double mp, bp, jp; /* corresponding probs */
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
               if (bq < conf->min_baseQ) {
                    base_skip = 1;
                    goto check_indel; /* FIXME: argh! */
               }
               if (toupper(ref_base) != 'N' && toupper(nt) != toupper(ref_base)) {
                    if (bq < conf->min_altbaseQ) {
                         base_skip = 1;
                         goto check_indel; /* FIXME: argh! */
                    }
               }
               /* FIXME test if this work as expected */
               /* FIXME should call cons here and test against it? */
               
               /* the following will correct base-pairs
                * down if they exceed the valid
                * sanger/phred limits. is it wise to do
                * this automatically? doesn't this
                * indicate a problem with the input ? */
               if (bq > 93) 
                    bq = 93; /* Sanger/Phred max */
               
               /* mapping quality. originally:
                * c = plp[i][j].b->core.qual + 33;
                * no need for check if mq is within user defined
                * limits. check was done in mplp_func */
               mq = p->b->core.qual;
               
               /* samtools check to detect Sanger max
                * value: problem is that an MQ Phred of
                * 255 means NA according to the samtools
                * spec (needed below). This however is not
                * detectable if the following original
                * samtools line gets executed, which is
                * why we remove it:
                * if (mq > 126) mq = 126;
                */
               
               /* "Merge" MQ and BQ if requested and if
                *  MAQP not 255 (not available):
                * P_jq = P_mq * + (1-P_mq) P_bq.
                */
               if ((conf->flag & MPLP_JOIN_BQ_AND_MQ) && mq != 255) {
                    /* Careful: q's are all ready to print
                     * sanger phred-scores. No need to do
                     * computation in phred-space as numbers
                     * won't get small enough.
                     */
                    mp = PHREDQUAL_TO_PROB(mq);
                    bp = PHREDQUAL_TO_PROB(bq);
                    /* precision?! */
                    jp = mp + (1.0 - mp) * bp;
                    jq = PROB_TO_PHREDQUAL(jp);
#ifdef DEBUG
                    LOG_DEBUG("P_M + (1-P_M) P_B:   %g + (1.0 - %g) * %g = %g  ==  Q%d + (1.0 - Q%d) * Q%d  =  Q%d\n",
                              mp+33, mp+33, bp+33, jp+33,  mq, mq, bq, jq);
#endif
                    final_q = jq;
               } else {
                    final_q = bq;
               }
               
               if (bam1_strand(p->b)) {
                    PLP_COL_ADD_BASEQUAL(& plp_col.qs_rv[nt4], final_q);
               } else {
                    PLP_COL_ADD_BASEQUAL(& plp_col.qs_fw[nt4], final_q);
               }
               
          } /* ! p->is_del */
          
          
     check_indel:
          
          /* coverage */
          if (p->is_del || p->is_refskip || 1 == base_skip) {
               num_skips += 1;
          }
          
          /* A pattern \+[0-9]+[ACGTNacgtn]+' indicates there is an
           * insertion between this reference position and the next
           * reference position. The length of the insertion is
           * given by the integer in the pattern, followed by the
           * inserted sequence. Similarly, a pattern
           * -[0-9]+[ACGTNacgtn]+ represents a deletion from the
           * reference. The deleted bases will be presented as ‘*’
           * in the following lines.
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
          
     }  /* end: for (j = 0; j < n_plp; ++j) { */
     
     /* FIXME: unused */
     f_depth = n_plp-num_skips;

     plp_col_mpileup_print(& plp_col, stdout);

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

    fprintf(stderr, "[%s] Note, the output differs from regular pileup (see http://samtools.sourceforge.net/pileup.shtml) in the following ways\n", __func__);
    fprintf(stderr, "[%s] - bases and qualities are merged into one field\n", __func__);
    fprintf(stderr, "[%s] - each base is immediately followed by its quality\n", __func__);
    fprintf(stderr, "[%s] - on request mapping and base call quality are merged (P_joined = P_mq + (1-P_mq)*P_bq\n", __func__);
    fprintf(stderr, "[%s] - indel events are removed from bases and qualities and summarized in an additional field\n", __func__);
    fprintf(stderr, "[%s] - reference matches are not replaced with , or .\n", __func__);

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


        process_plp(plp[i], n_plp[i], conf, ref, pos, ref_len, h->target_name[tid]);
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


void usage(const mplp_conf_t *mplp_conf) {
     fprintf(stderr, "Usage: %s [mpileup] [options] in.bam\n\n", MYNAME);
     fprintf(stderr, "Options:\n");
     /* regions */
     fprintf(stderr, "       -r STR       region in which pileup is generated [null]\n");
     fprintf(stderr, "       -l FILE      list of positions (chr pos) or regions (BED) [null]\n");
     /*  */
     fprintf(stderr, "       -d INT       max per-BAM depth to avoid excessive memory usage [%d]\n", mplp_conf->max_depth);
     fprintf(stderr, "       -f FILE      faidx indexed reference sequence file [null]\n");
     /* base call quality and baq */
     fprintf(stderr, "       -q INT       skip any base with baseQ smaller than INT [%d]\n", mplp_conf->min_baseQ);
     fprintf(stderr, "       -Q INT       skip nonref-bases with baseQ smaller than INT [%d]. Not active if ref is N\n", mplp_conf->min_altbaseQ);
     fprintf(stderr, "       -B           disable BAQ computation\n");
     fprintf(stderr, "       -E           extended BAQ for higher sensitivity but lower specificity\n");
     /* mapping quality */
     fprintf(stderr, "       -m INT       skip alignments with mapQ smaller than INT [%d]\n", mplp_conf->min_mq);
     fprintf(stderr, "       -M INT       cap mapping quality at INT [%d]\n", mplp_conf->max_mq);
     fprintf(stderr, "       -j           join mapQ and baseQ per base: P_e = P_mq + (1-P_mq) P_bq\n");
     /* misc */
     fprintf(stderr, "       -6           assume the quality is in the Illumina-1.3+ encoding\n");
     fprintf(stderr, "       -A           count anomalous read pairs\n");
}


int bam_mpileup(int argc, char *argv[])
{
    int c;
    int use_orphan = 0;
    mplp_conf_t mplp;

    memset(&mplp, 0, sizeof(mplp_conf_t));
    /* default pileup options */
    mplp.max_mq = 255; /* 60 */
    mplp.min_baseQ = 3; /* 13 */
    mplp.min_altbaseQ = 3; /* new */
    mplp.capQ_thres = 0;
    mplp.max_depth = 1000000; /* 250 */
    mplp.flag = MPLP_NO_ORPHAN | MPLP_REALN;
  
    while ((c = getopt(argc, argv, "r:l:d:f:Q:q:BEm:M:j6A:h")) >= 0) {
        switch (c) {
        case 'r': mplp.reg = strdup(optarg); break;
        case 'l': mplp.bed = bed_read(optarg); break;

        case 'd': mplp.max_depth = atoi(optarg); break;
        case 'f':
            mplp.fai = fai_load(optarg);
            if (mplp.fai == 0) 
                 return 1;
            break;

        case 'q': mplp.min_baseQ = atoi(optarg); break;
        case 'Q': mplp.min_altbaseQ = atoi(optarg); break;
        case 'B': mplp.flag &= ~MPLP_REALN; break;
        case 'E': mplp.flag |= MPLP_EXT_BAQ; break;

        case 'm': mplp.min_mq = atoi(optarg); break;
        case 'M': mplp.max_mq = atoi(optarg); break;
        case 'j': mplp.flag |= MPLP_JOIN_BQ_AND_MQ; break;

        case '6': mplp.flag |= MPLP_ILLUMINA13; break;
        case 'A': use_orphan = 1; break;

        case 'h': usage(& mplp); exit(0);
        case '?': fprintf(stderr, "ERROR: unrecognized arguments found. Exiting...\n"); exit(1);
        }
    }
    if (use_orphan) mplp.flag &= ~MPLP_NO_ORPHAN;
    
    assert(mplp.min_mq <= mplp.max_mq);
    assert(mplp.min_baseQ <= mplp.min_altbaseQ);
    

    if (argc == 1) {
        fprintf(stderr, "\n");
        usage(& mplp);
        return 1;
    }

    if (1 != argc - optind) {
        fprintf(stderr, "Need exactly one BAM file as last argument\n");
        return(EXIT_FAILURE);
    }

    mpileup(&mplp, 1, argv + optind);

    free(mplp.reg); free(mplp.pl_list);
    if (mplp.fai) fai_destroy(mplp.fai);
    if (mplp.bed) bed_destroy(mplp.bed);
    return 0;
}


int main(int argc, char **argv)
{
#ifdef BAM_NT_EXAMPLE
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


     if (argc>1 && strcmp(argv[1], "mpileup") == 0) {
          return bam_mpileup(argc-1, argv+1);
     } else {
          return bam_mpileup(argc, argv);
     }
}

