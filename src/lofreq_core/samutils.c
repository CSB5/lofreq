/* -*- c-file-style: "k&r"; indent-tabs-mode: nil; -*- */

/*********************************************************************
 *
 * FIXME update license
 *
 *********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <assert.h>

/* samtools includes */
#include "sam.h"
#include "kstring.h"

/* lofreq includes */
#include "log.h"
#include "vcf.h"
#include "plp.h"
#include "samutils.h"


char *
cigar_str_from_bam(const bam1_t *b)
{
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
/* cigar_str_from_bam() */


/* Count matches (OP_MATCH), mismatches (OP_MISMATCH), insertions
 * (OP_INS) and deletions (OP_DEL) for an aligned read. Written to
 * (preallocated, size 4) counts at indices given above. Returns non-0
 * on error. will ignore all non-del bases if their bq is below
 * min_bq.
 *
 * If quals is not NULL it will be used as a two dim array (has to be
 * preallocated) with OPs as first dim (len NUM_OP_CATS) and the
 * qualities of the bases as second dim. NOTE/FIXME: this uses bq for
 * mis/matches and INDEL_QUAL_DEFAULT for now in case of indels . The
 * number of elements corresponds to the count entry and can be at max
 * readlen.
 *
 * If target is non-NULL will ignore preloaded variant positions via
 * var_in_ign_list
 * 
 * Returns the total number of operations counted (excl clips or
 * bases<mq) or -1 on error
 */
int
count_cigar_ops(int *counts, int **quals,
                const bam1_t *b, const char *ref, int min_bq,
                char *target)
{
     const int INDEL_QUAL_DEFAULT = 40;
     int num_ops = 0;
     /* modelled after bam.c:bam_calend(), bam_format1_core() and
      * pysam's aligned_pairs (./pysam/csamtools.pyx)
      */
     uint32_t *cigar = bam1_cigar(b);
     const bam1_core_t *c = &b->core;
     uint32_t pos = c->pos; /* pos on genome */
     uint32_t qpos = 0; /* pos on read/query */
     uint32_t k, i;
#if 1
     int32_t qlen = (int32_t) bam_cigar2qlen(c, cigar); /* read length */
#else
     int qlen = b->core.l_qseq; /* read length */
#endif
     if (NULL==ref) {
          return -1;
     }

     assert(NULL != counts);
     memset(counts, 0, NUM_OP_CATS*sizeof(int));

     /* loop over cigar to get aligned bases
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
          if (op == BAM_CMATCH || op == BAM_CDIFF) {
               for (i=pos; i<pos+l; i++) {                             
                    int actual_op;
                    assert(qpos < qlen);
                    /* case agnostic */
                    char ref_nt = toupper(ref[i]);
                    char read_nt = bam_nt16_rev_table[bam1_seqi(bam1_seq(b), qpos)];
                    int bq = bam1_qual(b)[qpos];
#if 0
                    printf("[M]MATCH qpos,i,ref,read = %d,%d,%c,%c\n", qpos, i, ref_nt, read_nt);
#endif                    

                    if (bq<min_bq) {
                         /* fprintf(stderr, "M: ignoring because of bq=%d at %d (qpos %d)\n", bq, pos, qpos); */
                         qpos += 1;
                         continue;
                    }
#ifdef USE_SOURCEQUAL
                    if (target) {
                         var_t fake_var;
                         memset(&fake_var, 0, sizeof(var_t))
                         fake_var.chrom = target;
                         fake_var.pos = i;
                         /* FIXME evil hack. only work because var_in_ign_list only uses chrom and pos */
                         if (var_in_ign_list(&fake_var)) {
                              qpos += 1;
                              continue;
                         } 
                    }
#endif 
                    if (ref_nt != read_nt || op == BAM_CDIFF) {
                         actual_op = OP_MISMATCH;
                    } else {
                         actual_op = OP_MATCH;
                    }

 
                    counts[actual_op] += 1;
                    if (quals) {
                         quals[actual_op][counts[actual_op]-1] = bq;
                    }

                    qpos += 1;
               }
               pos += l;

          } else if (op == BAM_CINS) {
               for (i=pos; i<pos+l; i++) {
                    assert(qpos < qlen);
                    int bq = bam1_qual(b)[qpos];
                    if (bq<min_bq) {
                         /* fprintf(stderr, "I: ignoring because of bq=%d at %d (qpos %d)\n", bq, pos, qpos); */
                         qpos += 1;
                         continue;
                    }
#ifdef USE_SOURCEQUAL
                    if (target) {
                         var_t fake_var;
                         fake_var.chrom = target;
                         fake_var.pos = i; /* FIXME does i make sense here, i.e. matches dbSNP entry? */
                         if (var_in_ign_list(&fake_var)) {
                              qpos += 1;
                              continue;
                         }
                    }
#endif

#if 0
                    printf("INS qpos,i = %d,None\n", qpos);
#endif
                    counts[OP_INS] += 1;
                    if (quals) {
                         quals[OP_INS][counts[OP_INS]-1] = INDEL_QUAL_DEFAULT; /* FIXME use iq */
                    }

                    qpos += 1;
               }
               
          } else if (op == BAM_CDEL || op == BAM_CREF_SKIP) {
               for (i=pos; i<pos+l; i++) {
#if 0
                    printf("DEL qpos,i = None,%d\n", i);
#endif

                    if (op == BAM_CDEL) {
#ifdef USE_SOURCEQUAL
                         if (target) {
                              var_t fake_var;
                              fake_var.chrom = target;
                              fake_var.pos = i; /* FIXME does i that make sense here, i.e. matches dbSNP entry? */
                              if (var_in_ign_list(&fake_var)) {
                                   continue;
                              }
                         }
#endif
                        counts[OP_DEL] += 1;
                        if (quals) {
                             quals[OP_DEL][counts[OP_DEL]-1] = INDEL_QUAL_DEFAULT; /* FIXME use dq */
                        }
                    }
               }
               pos += l;
               /* deletion: don't increase qpos */

          } else if (op == BAM_CSOFT_CLIP) {
#if 0
               printf("SOFT CLIP qpos = %d\n", qpos);
#endif
               qpos += l;

          } else if (op != BAM_CHARD_CLIP) {
               LOG_WARN("Unknown op %d in cigar %s\n", op, cigar_str_from_bam(b));

          }
     } /* for k */
     assert(pos == bam_calend(&b->core, bam1_cigar(b))); /* FIXME correct assert? what if hard clipped? */
     if (qpos != qlen) {
               LOG_FIXME("got qpos=%d and qlen=%d for cigar %s l_qseq %d\n", qpos, qlen, cigar_str_from_bam(b), b->core.l_qseq);
     }
     assert(qpos == qlen); /* FIXME correct assert? What if hard clipped? */

     num_ops = 0;
     for (i=0; i<NUM_OP_CATS; i++) {
          num_ops += counts[i];
     }
     return num_ops;
}
/* count_cigar_ops() */

