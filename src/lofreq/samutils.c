/* -*- c-file-style: "k&r"; indent-tabs-mode: nil; -*- */
/*********************************************************************
* The MIT License (MIT)
* 
* Copyright (c) 2013,2014 Genome Institute of Singapore
* 
* Permission is hereby granted, free of charge, to any person
* obtaining a copy of this software and associated documentation files
* (the "Software"), to deal in the Software without restriction,
* including without limitation the rights to use, copy, modify, merge,
* publish, distribute, sublicense, and/or sell copies of the Software,
* and to permit persons to whom the Software is furnished to do so,
* subject to the following conditions:
* 
* The above copyright notice and this permission notice shall be
* included in all copies or substantial portions of the Software.
* 
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
* MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
* NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
* BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
* ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
* CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
* SOFTWARE.
*
************************************************************************/


#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <assert.h>

/* htslib includes */
#include "htslib/sam.h"
#include "htslib/kstring.h"

/* lofreq includes */
#include "log.h"
#include "vcf.h"
#include "plp.h"
#include "samutils.h"

#if 0
/* libbam:bamaux.c */
extern void bam_init_header_hash(bam_hdr_t *header);
extern void bam_destroy_header_hash(bam_hdr_t *header);
#endif

#define INDEL_QUAL_DEFAULT 45

#define BUF_SIZE 1024

#define MAX_READ_LEN 8192

#ifdef USE_ALNERRPROF

void
free_alnerrprof(alnerrprof_t *alnerrprof)
{
     int i;
     for (i=0; i<alnerrprof->num_targets; i++) {
          if (alnerrprof->prop_len[i]){
               free(alnerrprof->props[i]);/* free right-away if no data */
          }
     }
     free(alnerrprof->props);
     free(alnerrprof->prop_len);
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
parse_alnerrprof_statsfile(alnerrprof_t *alnerrprof, const char *path, bam_hdr_t *bam_header)
{
     char line[BUF_SIZE];
     int i;
     int *max_obs_pos;
     const int default_read_len = 250;
     int free_bam_header_hash = 0;
     int rc;
     FILE *in = fopen(path, "r");

#if 0
     // HTSlib does not have ->hash or bam_init_header_hash().
     // Most likely all this is no longer needed.
     /* needed for finding tid from tname */
     if (bam_header->hash == 0) {
          bam_init_header_hash(bam_header);             
          free_bam_header_hash = 1;
     }
#endif

     max_obs_pos = calloc(bam_header->n_targets, sizeof(int));
     
     alnerrprof->num_targets = bam_header->n_targets;
     alnerrprof->prop_len = calloc(alnerrprof->num_targets, sizeof(int));
     alnerrprof->props = calloc(alnerrprof->num_targets, sizeof(double *));     
     for (i=0; i<alnerrprof->num_targets; i++) {
          alnerrprof->prop_len[i] = default_read_len;/* default alloc here and realloc later */
          alnerrprof->props[i] = calloc(alnerrprof->prop_len[i], sizeof(double));
     }
     i=-1; /* make sure value is not reused by accident; triggers clang warning though */

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

         tid = bam_name2id(bam_header, tname);
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
               free(alnerrprof->props[i]);/* no data for this tid: free */
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

#if 0
     free_bam_header_hash = 0; /* FIXME segfaults often for unknown reason */
     if (free_bam_header_hash) {
          bam_destroy_header_hash(bam_header);
     }
#endif
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
     uint32_t *cigar = bam_get_cigar(b);
     uint32_t k, i;
     const bam1_core_t *c = &b->core;
#if 0
     int32_t qlen = (int32_t) bam_cigar2qlen(c, cigar); /* read length */
#else
     int qlen = b->core.l_qseq; /* read length */
#endif
     uint32_t pos = c->pos; /* pos on genome */
     uint32_t qpos = 0; /* pos on read/query */
     uint32_t qpos_org = bam_is_rev(b) ? qlen-qpos-1 : qpos;/* original qpos before mapping as possible reverse */


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
                    char ref_nt = ref[i];
                    char read_nt = seq_nt16_str[bam_seqi(bam_get_seq(b), qpos)];
                    int bq = bam_get_qual(b)[qpos];
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
                    qpos_org = bam_is_rev(b) ? qlen-qpos-1 : qpos;
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
                    qpos_org = bam_is_rev(b) ? qlen-qpos-1 : qpos;
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
               qpos_org = bam_is_rev(b) ? qlen-qpos-1 : qpos;

          } else if (op != BAM_CHARD_CLIP) {
               LOG_WARN("Unknown op %d in cigar %s\n", op, cigar_str_from_bam(b));

          }
     } /* for k */
     assert(pos == bam_calend(&b->core, bam_get_cigar(b))); /* FIXME correct assert? what if hard clipped? */
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



/* from char *bam_format1_core(const bam_hdr_t *header, const
 * bam1_t *b, int of) 
 */
char *
cigar_str_from_bam(const bam1_t *b)
{
     const bam1_core_t *c = &b->core;
     kstring_t str;
     int i;
     str.l = str.m = 0; str.s = 0;
     for (i = 0; i < c->n_cigar; ++i) {
          kputw(bam_get_cigar(b)[i]>>BAM_CIGAR_SHIFT, &str);
          kputc("MIDNSHP=X"[bam_get_cigar(b)[i]&BAM_CIGAR_MASK], &str);
     }
     return str.s;
}
/* cigar_str_from_bam() */



/* Count matches (OP_MATCH), mismatches (OP_MISMATCH), insertions
 * (OP_INS) and deletions (OP_DEL) for an aligned read. Written to
 * (preallocated, size 4) counts at indices given above. Will ignore
 * all mis-/match bases if their bq is below min_bq.
 *
 * Returns the total number of operations counted (excl. clipped bases
 * or those with bq<min_bq) or -1 on error. Consecutive indels are
 * counted as one operation, using INDEL_QUAL_DEFAULT, which is
 * suboptimal. 0 is a valid return value, e.g. if all bases are below
 * the quality threshold.
 *
 * If quals is not NULL it will be used as a two dim array (has to be
 * preallocated) with OPs as first dim (len NUM_OP_CATS) and the
 * qualities of the bases as second dim. NOTE/FIXME: this uses bq for
 * mis/matches and INDEL_QUAL_DEFAULT for now in case of indels. The
 * number of elements corresponds to the count entry and can be at max
 * readlen.
 * 
 * If target is non-NULL will ignore preloaded variant positions via
 * var_in_ign_list
 *
 * WARNING code duplication with calc_read_alnerrprof but merging the
 * two functions was too complicated (and the latter is unused anyway)
 */
int
count_cigar_ops(int *counts, int **quals, const bam1_t *b,
                const char *ref, int min_bq, char *target)
{
#if 0
#define TRACE 1
#endif
     int num_ops = 0;
     /* modelled after bam.c:bam_calend(), bam_format1_core() and
      * pysam's aligned_pairs (./pysam/csamtools.pyx)
      */
     uint32_t *cigar = bam_get_cigar(b);
     const bam1_core_t *c = &b->core;
     uint32_t tpos = c->pos; /* pos on genome */
     uint32_t qpos = 0; /* pos on read/query */
     uint32_t k, i;
#if 0
     int32_t qlen = (int32_t) bam_cigar2qlen(c, cigar); /* read length */
#else
     int qlen = b->core.l_qseq; /* read length */
#endif

     if (! ref) {
          return -1;
     }
     if (! counts) {
          return -1;
     }

     memset(counts, 0, NUM_OP_CATS*sizeof(int));

     /* loop over cigar to get aligned bases
      *
      * read: bam_format1_core(NULL, b, BAM_OFDEC);
      */
     for (k=0; k < c->n_cigar; ++k) { /* n_cigar: number of cigar operations */
          int op = cigar[k] & BAM_CIGAR_MASK; /* the cigar operation */
          uint32_t l = cigar[k] >> BAM_CIGAR_SHIFT;

          /* following conditionals could be collapsed to much shorter
           * code, but we keep them roughly as they were in pysam's
           * aligned_pairs to make later comparison and handling of
           * indels easier
           */
          if (op == BAM_CMATCH || op == BAM_CDIFF) {
               for (i=tpos; i<tpos+l; i++) {                             
                    int actual_op;
                    assert(qpos < qlen);
                    char ref_nt = ref[i];
                    char read_nt = seq_nt16_str[bam_seqi(bam_get_seq(b), qpos)];
                    int bq = bam_get_qual(b)[qpos];

                    if (ref_nt != read_nt || op == BAM_CDIFF) {
                         actual_op = OP_MISMATCH;
                    } else {
                         actual_op = OP_MATCH;
                    }

                    /* ignoring base if below min_bq, independent of type */
                    if (bq<min_bq) {
#ifdef TRACE
                         fprintf(stderr, "TRACE(%s): [M]MATCH ignoring base because of bq=%d at %d (qpos %d)\n", bam_get_qname(b), bq, i, qpos);
#endif
                         qpos += 1;
                         continue;
                    }

                    /* for mismatches only */
                    if (target && actual_op == OP_MISMATCH) {
                         var_t fake_var;
                         memset(&fake_var, 0, sizeof(var_t));
                         fake_var.chrom = target;
                         fake_var.pos = i;
                         /* FIXME evil, evil hack. only works as long as var_in_ign_list only uses chrom and pos */
                         if (var_in_ign_list(&fake_var)) {

#ifdef TRACE
                              fprintf(stderr, "TRACE(%s): MM: ignoring because in ign list at %d (qpos %d)\n", bam_get_qname(b), i, qpos);
#endif
                              qpos += 1;
                              continue;
                         } 
                    }

#ifdef TRACE
                    fprintf(stderr, "TRACE(%s): adding [M]MATCH qpos,tpos,ref,read,bq = %d,%d,%c,%c,%d\n", bam_get_qname(b), qpos, tpos, ref_nt, read_nt, bq);
#endif                    
                    counts[actual_op] += 1;
                    if (quals) {
                         quals[actual_op][counts[actual_op]-1] = bq;
                    }

                    qpos += 1;
               }
               tpos += l;

          } else if (op == BAM_CINS || op == BAM_CDEL) {

               if (target) {
                    /* vcf: 
                     * indel at tpos 1 means, that qpos 2 is an insertion  (e.g. A to AT)
                     * del at tpos 1 means, that qpos 2 is missing (e.g. AT to A)
                     */
                    var_t fake_var;
                    fake_var.chrom = target;
                    fake_var.pos = tpos;
                    if (op==BAM_CINS) {
                         fake_var.pos -= 1;
                    }
                    /* FIXME see above: only works as long as var_in_ign_list only uses chrom and pos */
                    if (var_in_ign_list(&fake_var)) {
                         if (op == BAM_CINS) {
                              qpos += l;
                         }
#ifdef TRACE
                         fprintf(stderr, "TRACE(%s): %c: ignoring because in ign list at tpos %d (qpos %d)\n", bam_get_qname(b), op == BAM_CINS? 'I':'D', tpos, qpos);
#endif
                         continue;
                    }
               }

#ifdef TRACE
               fprintf(stderr, "TRACE(%s): adding %c qpos,tpos = %d,%d\n", bam_get_qname(b), op==BAM_CINS?'I':'D', qpos, tpos);
#endif                    

               if (op == BAM_CINS) {
                    counts[OP_INS] += 1; /* counts indel as 1 operation only */
                    if (quals) {
                         quals[OP_INS][counts[OP_INS]-1] = INDEL_QUAL_DEFAULT; /* FIXME use iq */
                    }
                    qpos += l;/* forward query pos by length of operation */

               } else if (op == BAM_CDEL) {
                    counts[OP_DEL] += 1; /* counts indel as 1 operation only */
                    if (quals) {
                         quals[OP_DEL][counts[OP_DEL]-1] = INDEL_QUAL_DEFAULT; /* FIXME use dq */
                    }
                    tpos += l; /* forward genome pos by length of operation */

               } else {
                    LOG_FATAL("%s\n", "INTERNAL ERROR: should never get here");
                    exit(1);
               }

          } else if (op == BAM_CREF_SKIP) {
               tpos += l;

          } else if (op == BAM_CSOFT_CLIP) {
#if 0
               printf("SOFT CLIP qpos = %d\n", qpos);
#endif
               qpos += l;

          } else if (op != BAM_CHARD_CLIP) {
               LOG_WARN("Untested op %d in cigar %s\n", op, cigar_str_from_bam(b));
               /* don't think we need to do anything here */
          }
     } /* for k */

     assert(qpos == bam_calend(&b->core, bam_get_cigar(b))); /* FIXME correct assert? what if hard clipped? */
     if (qpos != qlen) {
          LOG_WARN("got qpos=%d and qlen=%d for cigar %s l_qseq %d in read %s\n", qpos, qlen, cigar_str_from_bam(b), b->core.l_qseq, bam_get_qname(b));
     }
     assert(qpos == qlen);

     num_ops = 0;
     for (i=0; i<NUM_OP_CATS; i++) {
          num_ops += counts[i];
#ifdef TRACE
          int j;
          for (j=0; j<counts[i]; j++) {
               fprintf(stderr, "TRACE(%s) op %s #%d: %d\n", bam_get_qname(b), op_cat_str[i], j, quals[i][j]);
          }
#endif
     }
     return num_ops;
}
/* count_cigar_ops() */
#undef TRACE


/* check match between reference and bam files. prints an error
 * message and return non-zero on mismatch 
*/
int checkref(char *fasta_file, char *bam_file)
{
     int i = -1;
     bam_hdr_t *header;
     faidx_t *fai;
     char *ref;
     int ref_len = -1;
     samFile *bam_fp;
     
     if (! file_exists(fasta_file)) {
          LOG_FATAL("Fsata file %s does not exist. Exiting...\n", fasta_file);
          return 1;
     }     

     if (0 != strcmp(bam_file, "-")  && ! file_exists(bam_file)) {
          LOG_FATAL("BAM file %s does not exist. Exiting...\n", bam_file);
          return 1;
     }     

     bam_fp = sam_open(bam_file, "r");
     header = sam_hdr_read(bam_fp);
     if (!header) {
          LOG_FATAL("Failed to read BAM header from %s\n", bam_file);
          return 1;
     }
     
     fai = fai_load(fasta_file);
     if (!fai) {
          LOG_FATAL("Failed to fasta index for %s\n", fasta_file);
          return 1;
     }
     
     for (i=0; i < header->n_targets; i++) {
          LOG_DEBUG("BAM header target %d of %d: name=%s len=%d\n", 
                    i+1, header->n_targets, header->target_name[i], header->target_len[i]);
          
          ref = faidx_fetch_seq(fai, header->target_name[i], 
                                0, 0x7fffffff, &ref_len);
          if (NULL == ref) {
               LOG_FATAL("Failed to fetch sequence %s from fasta file\n", header->target_name[i]);
               return -1;
          }
          if (header->target_len[i] != ref_len) {
               LOG_FATAL("Sequence length mismatch for sequence %s (%dbp in fasta; %dbp in bam)\n", 
                         header->target_name[i], header->target_len[i], ref_len);
               return -1;
          }
          free(ref);
     }
     
     fai_destroy(fai);
     bam_hdr_destroy(header);
     sam_close(bam_fp);

     return 0;
}
