/* -*- c-file-style: "k&r" -*-
 * http://emacswiki.org/emacs/IndentingC
 * http://www.emacswiki.org/emacs/LocalVariables
 * http://en.wikipedia.org/wiki/Indent_style 
 *
 * This is a modified version of samtools mpileup. The source code is
 * based on bam_plcmd_mod.c.
 *
 * The original source can be compiled with -DORIG
 * The main function in the new source is activated with -DMAIN
 */

#include <math.h>
#include <stdio.h>
#include <unistd.h>
#include <ctype.h>
#include <string.h>
#include <errno.h>
#include <assert.h>

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
    int min_refbaseQ, min_altbaseQ; /* new */
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


const char *bam_nt4_rev_table = "ACGTN"; /* as bam_nt16_rev_table */
#define NUM_NT4 5 /* strlen(bam_nt4_rev_table); */

#define MYNAME "lofreq_mpileup"
#define PHREDQUAL_TO_PROB(phred) (pow(10.0, -1.0*(phred)/10.0))
#define PROB_TO_PHREDQUAL(prob) ((int)(-10.0 * log10(prob)))


typedef struct {
     unsigned long int n; /* number of elements stored */
     int *data; /* actual array of data */

     size_t grow_by_size; /* if needed grow array by this value. will double previous size if 0 */
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

void int_varray_add_value(int_varray_t *a, int value)
{
    assert(NULL != a);
    if (a->n * sizeof(int) == a->alloced) {
        size_t size_to_alloc;
        
        if (0 == a->grow_by_size) {
             assert(SIZE_MAX - a->alloced > a->alloced);
             size_to_alloc = 0==a->n ? sizeof(int) : a->alloced*2;
        } else {
             assert(SIZE_MAX - a->alloced > a->grow_by_size);
             size_to_alloc = a->alloced + a->grow_by_size;
        }
#if 0
        fprintf(stderr, "DEBUG(%s:%s:%d): size=%lu alloced=%zd size_to_alloc=%zd sizeof(data)=%zd\n",
                __FILE__, __FUNCTION__, __LINE__, 
                a->n, a->alloced, size_to_alloc, sizeof(a->data));
#endif
        a->data = realloc(a->data, size_to_alloc);
        a->alloced = size_to_alloc;
    }

    a->data[a->n] = value;
    a->n++;
}


typedef struct {
     char *target;
     int pos;
     char ref_base;
     int_varray_t qs_fw[NUM_NT4];
     int_varray_t qs_rv[NUM_NT4];
     int heads;
     int tails;
} plp_col_t;


void plp_col_init(plp_col_t *p) {
    int i;

    p->target =  NULL;
    p->pos = -1;
    p->ref_base = '\0';
    for (i=0; i<NUM_NT4; i++) {
         int_varray_init(& p->qs_fw[i], 0);
         int_varray_init(& p->qs_rv[i], 0);
    }
    p->heads = 0;
    p->tails = 0;
}

void plp_col_free(plp_col_t *p) {
    int i;

    free(p->target);
    for (i=0; i<NUM_NT4; i++) {
         int_varray_free(& p->qs_fw[i]);
         int_varray_free(& p->qs_rv[i]);
    }
}

void plp_col_print(const plp_col_t *p, FILE *stream) {
     int i;
     
     fprintf(stream, "%s\t%d\t%c base-counts (fw/rv): ", 
             p->target, p->pos+1, p->ref_base);
     for (i=0; i<NUM_NT4; i++) {
          fprintf(stream, " %c:%lu/%lu",
                  bam_nt4_rev_table[i],
                  p->qs_fw[i].n,
                  p->qs_rv[i].n);
     }
     fprintf(stream, " heads:%d tails:%d\n", p->heads, p->tails);
     fprintf(stream, "\n");
}


/* FIXME also in lofreq_core:utils.c */
int file_exists(char *fname) 
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

static inline void pileup_seq(const bam_pileup1_t *p, int pos, int ref_len, const char *ref)
{
	int j;
	if (p->is_head) {
		putchar('^');
		putchar(p->b->core.qual > 93? 126 : p->b->core.qual + 33);
	}
	if (!p->is_del) {
		int c = bam_nt16_rev_table[bam1_seqi(bam1_seq(p->b), p->qpos)];
		if (ref) {
			int rb = pos < ref_len? ref[pos] : 'N';
			if (c == '=' || bam_nt16_table[c] == bam_nt16_table[rb]) c = bam1_strand(p->b)? ',' : '.';
			else c = bam1_strand(p->b)? tolower(c) : toupper(c);
		} else {
			if (c == '=') c = bam1_strand(p->b)? ',' : '.';
			else c = bam1_strand(p->b)? tolower(c) : toupper(c);
		}
		putchar(c);
	} else putchar(p->is_refskip? (bam1_strand(p->b)? '<' : '>') : '*');
	if (p->indel > 0) {
		putchar('+'); printw(p->indel, stdout);
		for (j = 1; j <= p->indel; ++j) {
			int c = bam_nt16_rev_table[bam1_seqi(bam1_seq(p->b), p->qpos + j)];
			putchar(bam1_strand(p->b)? tolower(c) : toupper(c));
		}
	} else if (p->indel < 0) {
		printw(p->indel, stdout);
		for (j = 1; j <= -p->indel; ++j) {
			int c = (ref && (int)pos+j < ref_len)? ref[pos+j] : 'N';
			putchar(bam1_strand(p->b)? tolower(c) : toupper(c));
		}
	}
	if (p->is_tail) putchar('$');
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
		if (ret < 0) break;
		if (b->core.tid < 0 || (b->core.flag&BAM_FUNMAP)) { /* exclude unmapped reads */
			skip = 1;
			continue;
		}
		if (ma->conf->bed) { /* test overlap */
			skip = !bed_overlap(ma->conf->bed, ma->h->target_name[b->core.tid], b->core.pos, bam_calend(&b->core, bam1_cigar(b)));
			if (skip) continue;
		}
#ifdef ORIG
		if (ma->conf->rghash) { /* exclude read groups */
			uint8_t *rg = bam_aux_get(b, "RG");
			skip = (rg && bcf_str2id(ma->conf->rghash, (const char*)(rg+1)) >= 0);
			if (skip) continue;
		}
#endif
		if (ma->conf->flag & MPLP_ILLUMINA13) {
			int i;
			uint8_t *qual = bam1_qual(b);
			for (i = 0; i < b->core.l_qseq; ++i)
				qual[i] = qual[i] > 31? qual[i] - 31 : 0;
		}
		has_ref = (ma->ref && ma->ref_id == b->core.tid)? 1 : 0;
		skip = 0;
		if (has_ref && (ma->conf->flag&MPLP_REALN)) bam_prob_realn_core(b, ma->ref, (ma->conf->flag & MPLP_EXT_BAQ)? 3 : 1);
		if (has_ref && ma->conf->capQ_thres > 10) {
			int q = bam_cap_mapQ(b, ma->ref, ma->conf->capQ_thres);
			if (q < 0) skip = 1;
			else if (b->core.qual > q) b->core.qual = q;
		}
		else if (b->core.qual < ma->conf->min_mq) skip = 1; 
		else if ((ma->conf->flag&MPLP_NO_ORPHAN) && (b->core.flag&1) && !(b->core.flag&2)) skip = 1;
	} while (skip);
	return ret;
}


static int mpileup(mplp_conf_t *conf, int n, char **fn)
{
	mplp_aux_t **data;
	int i, tid, pos, *n_plp, tid0 = -1, beg0 = 0, end0 = 1u<<29, ref_len, ref_tid = -1, max_depth;
	const bam_pileup1_t **plp;
	bam_mplp_t iter;
	bam_header_t *h = 0;
	char *ref;
#ifdef ORIG
	void *rghash = 0;
#endif
	kstring_t buf;
	mplp_pileup_t gplp;

    /* create a fake struct, whose only member of interest here is an
     * int of value 1. originally bam_sample_t *sm. */
    typedef struct {
         int n;
    } bam_sample_dummy_t;
	bam_sample_dummy_t *sm;
	sm = calloc(1, sizeof(bam_sample_dummy_t));
    assert(NULL != sm);
	sm->n = 1; /* that's all we need sm for */
    

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

	memset(&gplp, 0, sizeof(mplp_pileup_t));
	memset(&buf, 0, sizeof(kstring_t));
#ifdef ORIG
	memset(&bc, 0, sizeof(bcf_call_t));
#endif
	data = calloc(n, sizeof(void*));
	plp = calloc(n, sizeof(void*));
	n_plp = calloc(n, sizeof(int*));

	fprintf(stderr, "[%s] Note, the output differs from regular pileup (see http://samtools.sourceforge.net/pileup.shtml) in the following ways\n", __func__);
	fprintf(stderr, "[%s] - bases and qualities are merged into one field\n", __func__);
    fprintf(stderr, "[%s] - each base is immediately followed by its quality\n", __func__);
	fprintf(stderr, "[%s] - on request mapping and base call quality are merged (P_joined = P_mq * + (1-P_mq) P_bq\n", __func__);
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
	gplp.n = sm->n;
	gplp.n_plp = calloc(sm->n, sizeof(int));
	gplp.m_plp = calloc(sm->n, sizeof(int));
	gplp.plp = calloc(sm->n, sizeof(void*));


	if (tid0 >= 0 && conf->fai) { /* region is set */
		ref = faidx_fetch_seq(conf->fai, h->target_name[tid0], 0, 0x7fffffff, &ref_len);
		ref_tid = tid0;
		for (i = 0; i < n; ++i) data[i]->ref = ref, data[i]->ref_id = tid0;
	} else ref_tid = -1, ref = 0;
	iter = bam_mplp_init(n, mplp_func, (void**)data);
	max_depth = conf->max_depth;
	if (max_depth * sm->n > 1<<20)
		fprintf(stderr, "(%s) Max depth is above 1M. Potential memory hog!\n", __func__);
	if (max_depth * sm->n < 8000) {
		max_depth = 8000 / sm->n;
		fprintf(stderr, "<%s> Set max per-file depth to %d\n", __func__, max_depth);
	}
	/* max_indel_depth = conf->max_indel_depth * sm->n; */
	bam_mplp_set_maxcnt(iter, max_depth);
	while (bam_mplp_auto(iter, &tid, &pos, n_plp, plp) > 0) {
        char ref_base;

		if (conf->reg && (pos < beg0 || pos >= end0)) continue; /* out of the region requested */
		if (conf->bed && tid >= 0 && !bed_overlap(conf->bed, h->target_name[tid], pos, pos+1)) continue;
		if (tid != ref_tid) {
			free(ref); ref = 0;
			if (conf->fai) ref = faidx_fetch_seq(conf->fai, h->target_name[tid], 0, 0x7fffffff, &ref_len);
			for (i = 0; i < n; ++i) data[i]->ref = ref, data[i]->ref_id = tid;
			ref_tid = tid;
		}

            ref_base = (ref && pos < ref_len)? ref[pos] : 'N';
     		printf("%s\t%d\t%c", h->target_name[tid], pos + 1, ref_base);

			for (i = 0; i < n; ++i) { /* NOTE n is fixed to 1 */
				int j;
                plp_col_t plp_col;
                plp_col_init(& plp_col);

                plp_col.target = strdup(h->target_name[tid]);
                plp_col.pos = pos;
                plp_col.ref_base = ref_base;

				printf("\t%d\t", n_plp[i]);
				if (n_plp[i] == 0) {
                     printf("*\t*"); /* FIXME: printf() is very slow... */
				} else {
#ifdef ORIG
					for (j = 0; j < n_plp[i]; ++j)
						pileup_seq(plp[i] + j, pos, ref_len, ref);
					putchar('\t');
					for (j = 0; j < n_plp[i]; ++j) {
						const bam_pileup1_t *p = plp[i] + j;
						int c = bam1_qual(p->b)[p->qpos] + 33;
						if (c > 126) c = 126;
						putchar(c);
					}
					if (conf->flag & MPLP_PRINT_MAPQ) {
						putchar('\t');
						for (j = 0; j < n_plp[i]; ++j) {
							int c = plp[i][j].b->core.qual + 33;
							if (c > 126) c = 126;
							putchar(c);
						}
					}
				}
			}
			putchar('\n');
#else
                  
            /* FIXME separate parsing, fitlering and printing */
                    int num_heads = 0;
                    int num_tails = 0;
                    int num_ins = 0;
                    int sum_ins = 0;
                    int num_dels = 0;
                    int sum_dels = 0;

                    for (j = 0; j < n_plp[i]; ++j) {

                      /* merged modified pileup_seq() in here 
                       */
                      const bam_pileup1_t *p = plp[i] + j;
                      int nt, nt4;
                      int mq, bq, jq, final_q; /* sanger phred scores */
                      double mp, bp, jp; /* corresponding probs */

                      if (! p->is_del) {
                           if (p->is_head) {
                                plp_col.heads += 1;
                                num_heads += 1;
                           }
                           if (p->is_tail) {
                                plp_col.tails += 1;
                                num_tails += 1;
                           }

                           nt = bam_nt16_rev_table[bam1_seqi(bam1_seq(p->b), p->qpos)];
                           nt4 = bam_nt16_nt4_table[bam1_seqi(bam1_seq(p->b), p->qpos)];
                           nt = bam1_strand(p->b)? tolower(nt) : toupper(nt);
                           putchar(nt);
                           

                           /* FIXME get rid of unnecessary +- 33
                            * business here. Use only Phred values
                            * which are the scores -33.
                            */
                           bq = bam1_qual(p->b)[p->qpos] + 33;
                           if (bq > 126) bq = 126; /* Sanger max */

                           /* mapping quality. originally:
                            * conf->flag & MPLP_PRINT_MAPQ
                            * c = plp[i][j].b->core.qual + 33;
                            */
                           mq = p->b->core.qual + 33;
                           
                           /* samtools check to detect Sanger max
                            * value: problem is that an MQ Phred of
                            * 255 means NA according to the samtools
                            * spec (needed below). This however is not
                            * detectable if the following original
                            * samtools line gets executed:

                            if (mq > 126) mq = 126;
                           */

                           /* "Merge" MQ and BQ if requested and if
                            *  MAQP not 255 (not available):
                            * P_jq = P_mq * + (1-P_mq) P_bq.
                            */
                           if ((conf->flag & MPLP_JOIN_BQ_AND_MQ) && (mq-33) != 255) {
                                /* Careful: q's are all ready to print
                                 * sanger phred-scores. No need to do
                                 * computation in phred-space as numbers
                                 * won't get small enough.
                                 */
                                mp = PHREDQUAL_TO_PROB(mq - 33);
                                bp = PHREDQUAL_TO_PROB(bq - 33);
                                /* precision?! */
                                jp = mp + (1.0 - mp) * bp;
                                jq = PROB_TO_PHREDQUAL(jp) + 33;
#if DEBUG
                                printf("P_M + (1-P_M) P_B:   %g + (1.0 - %g) * %g = %g  ==  Q%d + (1.0 - Q%d) * Q%d  =  Q%d\n",
                                       mp, mp, bp, jp,  mq-33, mq-33, bq-33, jq-33);
#endif
                                putchar(jq);
                                final_q = jq;
                           } else {
                                final_q = bq;
                                putchar(bq);
                           }

                           if (bam1_strand(p->b)) {
                                int_varray_add_value(
                                     & plp_col.qs_rv[nt4], final_q);
                           } else {
                                int_varray_add_value(
                                     & plp_col.qs_fw[nt4], final_q);
                           }


                           /*putchar(' ');*/
                      } /* ! p->is_del */

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
                             /* originally: + */
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
                             num_ins += 1;
                             sum_ins += p->indel;
                          } else if (p->indel < 0) {
                             /* originally: - */
                             num_dels += 1;
                             sum_dels -= p->indel;
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
                          /* as an aside, if we wanted to include the
                             indels, do we need to keep track of heads
                             and tails as well?! */
                      }
                     
					} /* end: for (j = 0; j < n_plp[i]; ++j) { */
                    printf("\t#heads=%d #tails=%d #ins=%d ins_len=%.1f #del=%d del_len=%.1f\n",
                           num_heads, num_tails,
                           num_ins, num_ins ? sum_ins/(float)num_ins : 0,
                           num_dels, num_dels ? sum_dels/(float)num_dels : 0);
#if PRINT_PLP_COL
                    plp_col_print(& plp_col, stderr);
#endif
                    plp_col_free(& plp_col);
				}
			}
#endif
		}
    free(sm); free(buf.s);
	for (i = 0; i < gplp.n; ++i) free(gplp.plp[i]);
	free(gplp.plp); free(gplp.n_plp); free(gplp.m_plp);

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


int bam_mpileup(int argc, char *argv[])
{
	int c;
    int use_orphan = 0;
	mplp_conf_t mplp;
	memset(&mplp, 0, sizeof(mplp_conf_t));
	mplp.max_mq = 255; /* 60 */
	mplp.min_refbaseQ = 3; /* min_baseQ = 13 */
	mplp.min_altbaseQ = 20; /* FIXME use */
    mplp.capQ_thres = 0;
	mplp.max_depth = 1000000; /* 250 */
	mplp.flag = MPLP_NO_ORPHAN | MPLP_REALN;
  
	while ((c = getopt(argc, argv, "Af:r:l:m:M:Bd:E6j")) >= 0) {
		switch (c) {
		case 'r': mplp.reg = strdup(optarg); break;
		case 'l': mplp.bed = bed_read(optarg); break;

		case 'd': mplp.max_depth = atoi(optarg); break;
		case 'f':
			mplp.fai = fai_load(optarg);
			if (mplp.fai == 0) return 1;
			break;

		case 'B': mplp.flag &= ~MPLP_REALN; break;
		case 'E': mplp.flag |= MPLP_EXT_BAQ; break;

		case 'm': mplp.min_mq = atoi(optarg); break;
		case 'M': mplp.max_mq = atoi(optarg); break;
		case 'j': mplp.flag |= MPLP_JOIN_BQ_AND_MQ; break;

		case '6': mplp.flag |= MPLP_ILLUMINA13; break;
		case 'A': use_orphan = 1; break;

        /*case 'Q': mplp.min_baseQ = atoi(optarg); break; */

        case '?': fprintf(stderr, "ERROR: unrecognized arguments found. Exiting...\n"); exit(1);
        }
	}
	if (use_orphan) mplp.flag &= ~MPLP_NO_ORPHAN;
	if (argc == 1) {
		fprintf(stderr, "\n");
		fprintf(stderr, "Usage: %s [mpileup] [options] in.bam\n\n", MYNAME);
		fprintf(stderr, "Options:\n");
        /* regions */
		fprintf(stderr, "       -r STR       region in which pileup is generated [null]\n");
		fprintf(stderr, "       -l FILE      list of positions (chr pos) or regions (BED) [null]\n");
        /*  */
		fprintf(stderr, "       -d INT       max per-BAM depth to avoid excessive memory usage [%d]\n", mplp.max_depth);
		fprintf(stderr, "       -f FILE      faidx indexed reference sequence file [null]\n");
        /* baq */
		fprintf(stderr, "       -B           disable BAQ computation\n");
		fprintf(stderr, "       -E           extended BAQ for higher sensitivity but lower specificity\n");
        /* mapping quality */
		fprintf(stderr, "       -m INT       skip alignments with mapQ smaller than INT [%d]\n", mplp.min_mq);
		fprintf(stderr, "       -M INT       cap mapping quality at INT [%d]\n", mplp.max_mq);
		fprintf(stderr, "       -C INT       parameter for adjusting mapQ; 0 to disable [0]\n");
        /* new options */
		fprintf(stderr, "       -j           join mapQ and baseQ per base: P_e = P_mq + (1-P_mq) P_bq\n");
        /* misc */
		fprintf(stderr, "       -6           assume the quality is in the Illumina-1.3+ encoding\n");
		fprintf(stderr, "       -A           count anomalous read pairs\n");
#ifdef ORIG
		fprintf(stderr, "       -R           ignore RG tags\n");
#endif
        /* orig -Q doesn't affect samtools mpileup */
		return 1;
	}

    if (1 != argc - optind) {
		fprintf(stderr, "Need exactly one BAM file as argument\n");
        return(EXIT_FAILURE);
    }
    mpileup(&mplp, 1, argv + optind);
#ifdef ORIG
	if (mplp.rghash) bcf_str2id_thorough_destroy(mplp.rghash);
#endif
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

     fprintf(stderr, "FIXME: make bq and mq quality filtering work\n");
#endif


     if (argc>1 && strcmp(argv[1], "mpileup") == 0) {
          return bam_mpileup(argc-1, argv+1);
     } else {
          return bam_mpileup(argc, argv);
     }
}

