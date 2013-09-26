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


/* Count matches (OP_MATCH), mismatches (OP_MISMATCH),
 * insertions (OP_INS) and deletions (OP_DEL) for an
 * aligned read. Written to (preallocated, size 4) counts at indices
 * given above. Returns non-0 on error. will ignore all non-del bases
 * if thei bq is below min_bq.
 *
 * A slightly more flexible version of this would return an int array
 * of qlen length where each position is set to the corresponding BAM
 * op
 */
int
count_matches(int *counts, const bam1_t *b, const char *ref, int min_bq)
{
     /* modelled after bam.c:bam_calend(), bam_format1_core() and
      * pysam's aligned_pairs (./pysam/csamtools.pyx)
      */
     uint32_t *cigar = bam1_cigar(b);
     const bam1_core_t *c = &b->core;
     uint32_t pos = c->pos; /* pos on genome */
     uint32_t qpos = 0; /* pos on read/query */
     uint32_t k, i;
#ifndef NDEBUG
     int32_t qlen = (int32_t) bam_cigar2qlen(c, bam1_cigar(b)); /* read length */
#endif
     if (NULL==ref) {
          return 1;
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
                    assert(qpos < qlen);
                    /* case agnostic */
                    char ref_nt = toupper(ref[i]);
                    char read_nt = bam_nt16_rev_table[bam1_seqi(bam1_seq(b), qpos)];
                    int bq = bam1_qual(b)[qpos];
                    if (bq<min_bq) {
                         /* fprintf(stderr, "M: ignoring because of bq=%d at %d (qpos %d)\n", bq, pos, qpos); */
                         qpos += 1;
                         continue;
                    }
#if 0
                    printf("[M]MATCH qpos,i,ref,read = %d,%d,%c,%c\n", qpos, i, ref_nt, read_nt);
#endif                    
                    if (ref_nt != read_nt || op == BAM_CDIFF) {
                         counts[OP_MISMATCH] += 1;
                    } else {
                         counts[OP_MATCH] += 1;
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
#if 0
                    printf("INS qpos,i = %d,None\n", qpos);
#endif
                    counts[OP_INS] += 1;
                    qpos += 1;
               }
               
          } else if (op == BAM_CDEL || op == BAM_CREF_SKIP) {
               for (i=pos; i<pos+l; i++) {
#if 0
                    printf("DEL qpos,i = None,%d\n", i);
#endif
                    if (op == BAM_CDEL) {
                        counts[OP_DEL] += 1;
                    }
               }
               pos += l;

          } else if (op == BAM_CSOFT_CLIP) {
               qpos += 1;

          } else if (op != BAM_CHARD_CLIP) {
               LOG_WARN("Unknow op %d in cigar %s\n", op, cigar_str_from_bam(b));
          }
     } /* for k */
     assert(pos == bam_calend(&b->core, bam1_cigar(b))); /* FIXME correct assert? what if hard clipped? */
     assert(qpos == qlen); /* FIXME correct assert? What if hard clipped? */

     return 0;
}
/* count_matches() */

