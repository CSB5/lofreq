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

/* libbam:bamaux.c */
extern void bam_init_header_hash(bam_header_t *header);
extern void bam_destroy_header_hash(bam_header_t *header);

#define INDEL_QUAL_DEFAULT 40

#define BUF_SIZE 1024



#ifdef USE_MAPERRPROF

void
free_alnerrprof(alnerrprof_t *alnerrprof)
{
     int i;
     free(alnerrprof->prop_len);
     for (i=0; i<alnerrprof->num_targets; i++) {
          if (alnerrprof->prop_len[i]){
               free(alnerrprof->props[i]);/* free right-away if no data */
          }
     }
     free(alnerrprof->props);
     alnerrprof->num_targets = -1;
}


void
normalize_alnerrprof(alnerrprof_t *alnerrprof)
{
     int i;

#if 0
     {/* fixme report */
          for (i=0; i<alnerrprof->num_targets; i++) {
               int j;
               fprintf(stderr, "FIXME in tid=%d len=%d: ", i, alnerrprof->prop_len[i]);
               for (j=0; j<alnerrprof->prop_len[i]; j++) {
                    fprintf(stderr, " %d:%g ", j, alnerrprof->props[i][j]);
               }
               fprintf(stderr, "\n");
          }
     }
#endif

     for (i=0; i<alnerrprof->num_targets; i++) {
          int j;
          double median = dbl_median(alnerrprof->props[i], alnerrprof->prop_len[i]);
#if 0
          fprintf(stderr, "FIXME tid=%d median=%g\n", i, median);
#endif
          for (j=0; j<alnerrprof->prop_len[i]; j++) {
               double val = alnerrprof->props[i][j] - median;
               if (val >= 0.0) {
                    alnerrprof->props[i][j]  = val;
               } else {
                    alnerrprof->props[i][j]  = 0.0;
               }

          }
     }
     
#if 0
     {/* fixme report */
          for (i=0; i<alnerrprof->num_targets; i++) {
               int j;
               fprintf(stderr, "FIXME out tid=%d len=%d: ", i, alnerrprof->prop_len[i]);
               for (j=0; j<alnerrprof->prop_len[i]; j++) {
                    fprintf(stderr, " %d:%g ", j, alnerrprof->props[i][j]);
               }
               fprintf(stderr, "\n");
          }
     }
#endif
}


/* will return non-0 on error. parsed error prof will be written to
 * alnerrprof. values are allocated here and should be freed with
 * free_alnerrprof */
int
parse_alnerrprof_statsfile(alnerrprof_t *alnerrprof, const char *path, bam_header_t *bam_header)
{
     char line[BUF_SIZE];
     int i;
     int *max_obs_pos;
     const int default_read_len = 250;
     int free_bam_header_hash = 0;
     int rc;
     FILE *in = fopen(path, "r");


     /* needed for finding tid from tname */
     if (bam_header->hash == 0) {
          bam_init_header_hash(bam_header);             
          free_bam_header_hash = 1;
     }

     max_obs_pos = calloc(bam_header->n_targets, sizeof(int));
     
     alnerrprof->num_targets = bam_header->n_targets;
     alnerrprof->prop_len = calloc(alnerrprof->num_targets, sizeof(int));
     alnerrprof->props = calloc(alnerrprof->num_targets, sizeof(double *));     
     for (i=0; i<alnerrprof->num_targets; i++) {
          alnerrprof->prop_len[i] = default_read_len;/* default alloc here and realloc later */
          alnerrprof->props[i] = calloc(alnerrprof->prop_len[i], sizeof(double));
     }
     i=-1; /* make sure i is not reused improperly; triggers clang warning though */

     while (NULL != fgets(line, BUF_SIZE, in)) {
          int pos = -1;
          char tname[BUF_SIZE];
          double prop = -1;
          unsigned long int count = -1;
          int tid = -1;
          if (line[0]=='#') {
               continue;
          }

          if (4 != sscanf(line, "%s\t%d\t%lg\t%lu\n", tname, &pos, &prop, &count)) {
              LOG_ERROR("Couldn't parse line %s\n", line);
              rc = 1;
              goto free_and_exit;
         }

         assert(prop>=0.0 && prop<=1.0);

         pos = pos - 1;
         assert(pos<MAX_READ_LEN);

         tid = bam_get_tid(bam_header, tname);
         if (-1 == tid) {
              LOG_ERROR("Target name '%s' found in error profile doesn't match any of the sequences in BAM header. Skipping and trying to continue...\n", tname);
              continue;
         }
         assert(tid<alnerrprof->num_targets);

         /* for later downsizing */
         if (pos+1 > max_obs_pos[tid]) {
              max_obs_pos[tid] = pos+1;
         }

         /* upsize if necessary */
         while (pos >= alnerrprof->prop_len[tid]) {
              LOG_DEBUG("upsizing pos+1=%d alnerrprof->prop_len[tid=%d]=%d\n\n", pos+1, tid, alnerrprof->prop_len[tid]);
              alnerrprof->prop_len[tid] *= 2;
              alnerrprof->props[tid] = realloc(alnerrprof->props[tid], alnerrprof->prop_len[tid] * sizeof(double));
         }
         alnerrprof->props[tid][pos] = prop;
     }

     /* downsize */
     for (i=0; i<alnerrprof->num_targets; i++) {
          if (max_obs_pos[i]) {
               LOG_DEBUG("downsizing alnerrprof->prop_len[tid=%d] to max %d\n", i, max_obs_pos[i]);
               alnerrprof->props[i] = realloc(alnerrprof->props[i], max_obs_pos[i] * sizeof(double));
          } else {
               free(alnerrprof->props[i]);/* free right-away if no data */
          }
          alnerrprof->prop_len[i] = max_obs_pos[i];
     }

#if 0
     {/* fixme report */
          for (i=0; i<alnerrprof->num_targets; i++) {
               int j;
               fprintf(stderr, "tid=%d len=%d: ", i, alnerrprof->prop_len[i]);
               for (j=0; j<alnerrprof->prop_len[i]; j++) {
                    fprintf(stderr, " %d:%g ", j, alnerrprof->props[i][j]);
               }
               fprintf(stderr, "\n");
               fprintf(stderr, "median for tid %d: %g for size %d\n",
                       i,
                       dbl_median(alnerrprof->props[i], alnerrprof->prop_len[i]),
                       alnerrprof->prop_len[i]);
          }
     }
#endif

     rc = 0;

free_and_exit:
     
     free(max_obs_pos);

     if (free_bam_header_hash) {
          bam_destroy_header_hash(bam_header);
     }
     fclose(in);

     return rc;
}



void
write_alnerrprof_stats(char *target_name, unsigned long int *alnerrprof_usedpos, 
                    double *alnerrprof, int max_obs_read_len, FILE *out)
{
     /* poor man's version (fw reads only and not taking quality into account):
      *  samtools view -h -F 0x10 $bam | samtools calmd -S -e - $reffa | cut -f 10 | awk '{for (i=0; i<length($0); i++) {if (substr($0, i, 1)!="=") {c[i]+=1}}} END {for (i in c) {print i, c[i], c[i]/NR}}' | sort -k 1 -n
      */
     int i;
     fprintf(out, "# Error, i.e. 'no-match' profile along read after subtracting base-call/indel quality\n");
     fprintf(out, "# Numbers are in scientific notation\n");
     fprintf(out, "# chrom\tread-pos\terror-freq\tcount\n");

     for (i=0; i<max_obs_read_len; i++) {
         double prop = 0.0;
         if (alnerrprof_usedpos[i]) {
              prop = alnerrprof[i]/(double)(alnerrprof_usedpos[i]);
         }
         fprintf(out, "%s\t%d\t%g\t%lu\n", target_name, i+1, prop, alnerrprof_usedpos[i]);/*, alnerrprof[i], alnerrprof_usedpos[i]);*/
     }
}



/* Counts probability of non-match count along the read after
 * subtracting error prob at that position (using the original
 * orientation). used_pos is an array of ints indicating whether
 * position was used or not (trimmed, clipped etc). alnerrprof and
 * used_pos must be of at least length b->core.l_qseq. Note: will add
 * to alnerrprof and used_pos, i.e. arrays should be initialized to 0 if
 * you don't want aggregate values.
 *
 * WARNING code duplication with count_cigar_ops but merging the two
 * functions is messy.
 */
void
calc_read_alnerrprof(double *alnerrprof, unsigned long int *used_pos, 
                   const bam1_t *b, const char *ref)
{
     /* modelled after bam.c:bam_calend(), bam_format1_core() and
      * pysam's aligned_pairs (./pysam/csamtools.pyx)
      */
     uint32_t *cigar = bam1_cigar(b);
     uint32_t k, i;
     const bam1_core_t *c = &b->core;
#if 0
     int32_t qlen = (int32_t) bam_cigar2qlen(c, cigar); /* read length */
#else
     int qlen = b->core.l_qseq; /* read length */
#endif
     uint32_t pos = c->pos; /* pos on genome */
     uint32_t qpos = 0; /* pos on read/query */
     uint32_t qpos_org = bam1_strand(b) ? qlen-qpos-1 : qpos;/* original qpos before mapping as possible reverse */


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
#if 0
                    printf("[M]MATCH qpos,i,ref,read = %d,%d,%c,%c\n", qpos, i, ref_nt, read_nt);
#endif                    

                    if (ref_nt != 'N') {
                         if (ref_nt != read_nt || op == BAM_CDIFF) {
                              alnerrprof[qpos_org] += (1.0 - PHREDQUAL_TO_PROB(bq));
                         } /* otherwise leave at 0.0 but count anyway */
                         used_pos[qpos_org] += 1;
                    }
                    qpos += 1;
                    qpos_org = bam1_strand(b) ? qlen-qpos-1 : qpos;
               }
               pos += l;

          } else if (op == BAM_CINS) {
               for (i=pos; i<pos+l; i++) {
                    assert(qpos < qlen);
                    
                    alnerrprof[qpos] += (1.0 - PHREDQUAL_TO_PROB(INDEL_QUAL_DEFAULT));
                    used_pos[qpos] += 1;
#if 0
                    printf("INS qpos,i = %d,None\n", qpos);
#endif
                    qpos += 1;
                    qpos_org = bam1_strand(b) ? qlen-qpos-1 : qpos;
               }
               
          } else if (op == BAM_CDEL || op == BAM_CREF_SKIP) {
               for (i=pos; i<pos+l; i++) {
#if 0
                    printf("DEL qpos,i = None,%d\n", i);
#endif

                    if (op == BAM_CDEL) {
                         alnerrprof[qpos] += (1.0 - PHREDQUAL_TO_PROB(INDEL_QUAL_DEFAULT));
                         used_pos[qpos] += 1;
                    }
               }
               pos += l;
               /* deletion: don't increase qpos */

          } else if (op == BAM_CSOFT_CLIP) {
#if 0
               printf("SOFT CLIP qpos = %d\n", qpos);
#endif
               qpos += l;
               qpos_org = bam1_strand(b) ? qlen-qpos-1 : qpos;

          } else if (op != BAM_CHARD_CLIP) {
               LOG_WARN("Unknown op %d in cigar %s\n", op, cigar_str_from_bam(b));

          }
     } /* for k */
     assert(pos == bam_calend(&b->core, bam1_cigar(b))); /* FIXME correct assert? what if hard clipped? */
     if (qpos != qlen) {
          LOG_FIXME("got qpos=%d and qlen=%d for cigar %s l_qseq %d\n", qpos, qlen, cigar_str_from_bam(b), b->core.l_qseq);
     }
     assert(qpos == qlen); /* FIXME correct assert? What if hard clipped? */

#if 0
     fprintf(stderr, "%s:", __FUNCTION__);
     for (i=0; i< b->core.l_qseq; i++) {
          fprintf(stderr, " %g/%d", alnerrprof[i], used_pos[i]);
     }
     fprintf(stderr, "\n");
#endif
}
#endif



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
 *
 * WARNING code duplication with count_cigar_ops but merging the two
 * functions was too complicated
 */
int
count_cigar_ops(int *counts, int **quals,
                const bam1_t *b, const char *ref, int min_bq,
                char *target)
{
     int num_ops = 0;
     /* modelled after bam.c:bam_calend(), bam_format1_core() and
      * pysam's aligned_pairs (./pysam/csamtools.pyx)
      */
     uint32_t *cigar = bam1_cigar(b);
     const bam1_core_t *c = &b->core;
     uint32_t pos = c->pos; /* pos on genome */
     uint32_t qpos = 0; /* pos on read/query */
     uint32_t k, i;
#if 0
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
                         memset(&fake_var, 0, sizeof(var_t));
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
                              fake_var.pos = i; /* FIXME does i make sense here, i.e. matches dbSNP entry? */
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

