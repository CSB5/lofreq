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

/* This file is partially based on samtools' bam_plcmd.c and very
 * likely needs an update whenever samtools/libbam is updated
 *
 */

#include <ctype.h>
#include <assert.h>
#include <errno.h>
#include <fenv.h>

#include "htslib/kstring.h"
#include "htslib/sam.h"

#include "log.h"
#include "plp.h"
#include "vcf.h"
#include "samutils.h"
#include "snpcaller.h"
#include "bam_md_ext.h"

const char *bam_nt4_rev_table = "ACGTN";


int missing_baq_warning_printed = 0;

/* from bedidx.c */
void *bed_read(const char *fn);
void bed_destroy(void *_h);
int bed_overlap(const void *_h, const char *chr, int beg, int end);

/* From the SAM spec: "tags starting with `X', `Y' and `Z' or tags
 * containing lowercase letters in either position are reserved for
 * local use".
*/
#define SRC_QUAL_TAG "sq"

/* results on icga dream syn1.2 suggest that somatic calls made extra
 * with this settings are likely fp whereas the ones missing a likely
 * tp, therefore disabled */
#undef MERGEQ_FOR_CONS_CALL


const unsigned char bam_nt4_table[256] = {
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

typedef struct {
     samFile *fp;
     hts_itr_t* iter;
     bam_hdr_t *h;
     int ref_id;
     char *ref;
     const mplp_conf_t *conf;
} mplp_aux_t;

typedef struct {
    int n;
    int *n_plp, *m_plp;
    bam_pileup1_t **plp;
} mplp_pileup_t;


#ifdef USE_ALNERRPROF
static alnerrprof_t *alnerrprof = NULL;
#endif


/* initialize members of preallocated mplp_conf */
void init_mplp_conf(mplp_conf_t *c)
{
     memset(c, 0, sizeof(mplp_conf_t));
     c->max_mq = DEFAULT_MAX_MQ;
     c->min_mq = DEFAULT_MIN_MQ;
     c->def_nm_q = DEFAULT_DEF_NM_QUAL;
     c->min_plp_bq = DEFAULT_MIN_PLP_BQ;/* note: different from DEFAULT_MIN_BQ */
     c->min_plp_idq = DEFAULT_MIN_PLP_IDQ;
     c->max_depth = DEFAULT_MAX_PLP_DEPTH;
     c->flag = MPLP_NO_ORPHAN | MPLP_BAQ | MPLP_EXT_BAQ | MPLP_IDAQ;
}



/* convenience function */
int
base_count(const plp_col_t *p, char base)
{
     int b = bam_nt4_table[(int)base];
     return p->fw_counts[b] + p->rv_counts[b];
}



void
plp_col_init(plp_col_t *p) {
    int i;

    const int grow_by_size = 16384;

    p->target =  NULL;
    p->pos = -INT_MAX;
    p->ref_base = '\0';
    p->cons_base[0] = 'N'; p->cons_base[1] = '\0';
    p->coverage_plp = 0;
    p->num_bases = 0;
    p->num_ign_indels = 0;
    p->num_non_indels = 0;
    for (i=0; i<NUM_NT4; i++) {
         int_varray_init(& p->base_quals[i], grow_by_size);
         int_varray_init(& p->baq_quals[i], grow_by_size);
         int_varray_init(& p->map_quals[i], grow_by_size);
         int_varray_init(& p->source_quals[i], grow_by_size);
#ifdef USE_ALNERRPROF
         int_varray_init(& p->alnerr_qual[i], grow_by_size);
#endif
         p->fw_counts[i] = 0;
         p->rv_counts[i] = 0;
    }

    p->num_heads = p->num_tails = 0;

    p->num_ins = p->sum_ins = 0;
    int_varray_init(& p->ins_quals, grow_by_size);
    int_varray_init(& p->ins_map_quals, grow_by_size);
    int_varray_init(& p->ins_source_quals, grow_by_size);
    p->ins_event_counts = NULL;

    p->num_dels = p->sum_dels = 0;
    int_varray_init(& p->del_quals, grow_by_size);
    int_varray_init(& p->del_map_quals, grow_by_size);
    int_varray_init(& p->del_source_quals, grow_by_size);
    p->del_event_counts = NULL;

    p->non_ins_fw_rv[0] = p->non_ins_fw_rv[1] = 0;
    p->non_del_fw_rv[0] = p->non_del_fw_rv[1] = 0;

    p->has_indel_aqs = 0;
    p->hrun = 0;
}


void
plp_col_free(plp_col_t *p) {
    int i;

    free(p->target);
    for (i=0; i<NUM_NT4; i++) {
         int_varray_free(& p->base_quals[i]);
         int_varray_free(& p->baq_quals[i]);
         int_varray_free(& p->map_quals[i]);
         int_varray_free(& p->source_quals[i]);
#ifdef USE_ALNERRPROF
         int_varray_free(& p->alnerr_qual[i]);
#endif
    }

    int_varray_free(& p->ins_quals);
    int_varray_free(& p->ins_map_quals);
    int_varray_free(& p->ins_source_quals);
    int_varray_free(& p->del_quals);
    int_varray_free(& p->del_map_quals);
    int_varray_free(& p->del_source_quals);

    destruct_ins_event_counts(&p->ins_event_counts);
    destruct_del_event_counts(&p->del_event_counts);
}


void plp_col_debug_print(const plp_col_t *p, FILE *stream)
{
     int i;

     fprintf(stream, "%s\t%d\t%c\t%s\tcounts:rv/fw",
             p->target, p->pos+1, p->ref_base, p->cons_base);
     for (i=0; i<NUM_NT4; i++) {
          fprintf(stream, " %c:%lu/%lu",
                  bam_nt4_rev_table[i],
                  p->fw_counts[i],
                  p->rv_counts[i]);
     }

     fprintf(stream, " heads:%d tails:%d", p->num_heads, p->num_tails);
     fprintf(stream, " ins:%d del:%d", p->num_ins, p->num_dels);
     fprintf(stream, " hrun=%d", p->hrun);
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
void
plp_col_mpileup_print(const plp_col_t *p, FILE *stream)
{
     int i, j;

     fprintf(stream, "%s\t%d\t%c\t%d\t",
             p->target, p->pos+1, p->ref_base,p->coverage_plp);
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
/* plp_col_mpileup_print() */



void
dump_mplp_conf(const mplp_conf_t *c, FILE *stream)
{
     fprintf(stream, "mplp options\n");
     fprintf(stream, "  max_mq       = %d\n", c->max_mq);
     fprintf(stream, "  min_mq       = %d\n", c->min_mq);
     fprintf(stream, "  flag         = %d\n", c->flag);

     fprintf(stream, "  flag & MPLP_NO_ORPHAN  = %d\n", c->flag & MPLP_NO_ORPHAN ? 1:0);
     fprintf(stream, "  flag & MPLP_BAQ        = %d\n", c->flag & MPLP_BAQ ? 1:0);
     fprintf(stream, "  flag & MPLP_REDO_BAQ   = %d\n", c->flag & MPLP_REDO_BAQ ? 1:0);
     fprintf(stream, "  flag & MPLP_EXT_BAQ    = %d\n", c->flag & MPLP_EXT_BAQ ? 1:0);
     fprintf(stream, "  flag & MPLP_IDAQ       = %d\n", c->flag & MPLP_IDAQ ? 1:0);
     fprintf(stream, "  flag & MPLP_REDO_IDAQ = %d\n", c->flag & MPLP_REDO_IDAQ ? 1:0);
     fprintf(stream, "  flag & MPLP_USE_SQ     = %d\n", c->flag & MPLP_USE_SQ ? 1:0);
     fprintf(stream, "  flag & MPLP_ILLUMINA13 = %d\n", c->flag & MPLP_ILLUMINA13 ? 1:0);

     fprintf(stream, "  max_depth    = %d\n", c->max_depth);
     fprintf(stream, "  min_plp_bq   = %d\n", c->min_plp_bq);
     fprintf(stream, "  min_plp_idq  = %d\n", c->min_plp_idq);
     fprintf(stream, "  def_nm_q     = %d\n", c->def_nm_q);
     fprintf(stream, "  reg          = %s\n", c->reg);
     fprintf(stream, "  fa           = %p\n", c->fa);
     /*fprintf(stream, "  fai          = %p\n", c->fai);*/
     fprintf(stream, "  bed          = %p\n", c->bed);
     fprintf(stream, "  cmdline      = %s\n", c->cmdline);
}
/* dump_mplp_conf() */





static var_hash_t *source_qual_ign_vars_hash = NULL; /* must be declared NULL ! */


int
var_in_ign_list(var_t *var) {
     char *key = NULL;
     var_hash_t *match = NULL;

     /* using key_pos_only i.e. chrom and pos only
      *
      * NOTE: source quality will pass down fake vars without ref and
      * alt so only vcf_var_key_pos_only will work!
      */
     vcf_var_key_pos_only(&key, var);
     HASH_FIND_STR(source_qual_ign_vars_hash, key, match);
     free(key);

     if (NULL == match) {
          return 0;
     } else {
          return 1;
     }
}


void
source_qual_free_ign_vars()
{
     var_hash_free_table(source_qual_ign_vars_hash);
}


/* FIXME ignore variants outside given region  (mplp_conf.reg)
 (on top of bed as well)
 * and use tabix API if indexed */
int
source_qual_load_ign_vcf(const char *vcf_path, void *bed)
{
     vcf_file_t vcf_file;
     const int read_only_passed = 0;
     unsigned int num_total_vars = 0;
     unsigned int num_kept_vars = 0;

     if (vcf_file_open(& vcf_file, vcf_path,
                      HAS_GZIP_EXT(vcf_path), 'r')) {
         LOG_ERROR("Couldn't open %s\n", vcf_path);
         return 1;
     }

     if (0 !=  vcf_skip_header(& vcf_file)) {
         LOG_WARN("%s\n", "vcf_skip_header() failed");
         return 1;
     }

    /* as in lofreq_vcfset.c
     * WARN: partial code duplication
     */
    while (1) {
         var_t *var;
         char *key;
         int rc;

         vcf_new_var(&var);
         rc = vcf_parse_var(& vcf_file, var);
         if (rc) {
              vcf_free_var(&var);
              break;
         }
         num_total_vars += 1;

         if (read_only_passed && ! VCF_VAR_PASSES(var)) {
              vcf_free_var(&var);
              continue;
         }

         if (bed && ! bed_overlap(bed, var->chrom, var->pos, var->pos+1)) {
              vcf_free_var(&var);
              continue;
         }

         /* using key_pos_only i.e. chrom and pos only */
         vcf_var_key_pos_only(&key, var);
         /* since we only need the key and no other info we do
          * not need to save the var (and save NULL instead) */
         vcf_free_var(&var);
         var_hash_add(& source_qual_ign_vars_hash, key, NULL);
         /* FIXME key used within vcf_free_var; dont free! */
    }

    num_kept_vars = HASH_COUNT(source_qual_ign_vars_hash);
    if (num_kept_vars) {
         LOG_VERBOSE("Ignoring %d variants for SQ computation after reading %s\n",
                     num_kept_vars, vcf_path);
    } else {
         LOG_WARN("None of the %d variants in %s were kept\n",
                  num_total_vars, vcf_path);
    }
    vcf_file_close(& vcf_file);

    return 0;
}



/* Estimate as to how likely it is that this read, given the mapping,
 * comes from this reference genome. P(r not from g|mapping) = 1 - P(r
 * from g).
 *
 * Use base-qualities and poisson-binomial dist, similar to core SNV
 * calling, but return prob instead of pvalue (and subtract one
 * mismatch which is the SNV we are checking for the benefit of doubt;
 * affectively also means that prob MM=1 eg prob MM=0). If
 * nonmatch_qual is < 0, then keep all qualities as they are, i.e.
 * don't replace mismatches with lower values. Rationale: higher SNV
 * quals, means higher chance SNVs are real, therefore higher prob.
 * read does not come from genome. Otherwise use this phredscore as
 * default.
 *
 * Assuming independence of errors is okay, because if they are not
 * independent, then the prediction made is conservative.
 *
 * If target is non-NULL will ignore SNPs via var_in_ign_list
 *
 * Returns -1 on error or if NA. otherwise source quality
 *
 */
int
source_qual(const bam1_t *b, const char *ref,
            const int nonmatch_qual, char *target, int min_bq)
{
     int op_counts[NUM_OP_CATS];
     int **op_quals = NULL;

     double *probvec = NULL;
     int num_non_matches = -1; /* non-matching operations */
     int orig_num_non_matches = -1;
     double *err_probs = NULL; /* error probs (qualities) passed down to snpcaller. one for each op, no matter if matching or not */
     int num_err_probs; /* #elements in err_probs */

     long double unused_pval;
     int src_qual = -1;
     double src_prob = -1; /* prob of this read coming from genome */
     int err_prob_idx;
     int i, j;

     int qlen = b->core.l_qseq;

     /* alloc op_quals
      */
     if (NULL == (op_quals = malloc(NUM_OP_CATS * sizeof(int *)))) {
          fprintf(stderr, "FATAL: couldn't allocate memory at %s:%s():%d\n",
                  __FILE__, __FUNCTION__, __LINE__);
          exit(1);
     }
     for (i=0; i<NUM_OP_CATS; i++) {
          if (NULL == (op_quals[i] = malloc(qlen * sizeof(int)))) {/* over allocating */
               fprintf(stderr, "FATAL: couldn't allocate memory at %s:%s():%d\n",
                       __FILE__, __FUNCTION__, __LINE__);
               free(op_quals);
               exit(1);
          }
     }

     /* count match operations and get qualities for them
      */
     /* LOG_FIXME("%s\n", "Don't know ref name in count_cigar_ops which would be needed as hash key");*/
     num_err_probs = count_cigar_ops(op_counts, op_quals, b, ref, min_bq, target);
     if (1 > num_err_probs) {
#ifdef TRACE
          LOG_DEBUG("count_cigar_ops returned %d counts on read %s\n", num_err_probs, bam_get_qname(b));
#endif
          src_qual = -1;
          goto free_and_exit;
     }

     /* alloc and fill err_probs with quals returned per op-cat from
      * count_cigar_ops
      */
     if (NULL == (err_probs = malloc(num_err_probs * sizeof(double)))) {
          fprintf(stderr, "FATAL: couldn't allocate memory at %s:%s():%d\n",
                  __FILE__, __FUNCTION__, __LINE__);
          exit(1);
     }
     num_non_matches = 0;
     err_prob_idx = 0;
     for (i=0; i<NUM_OP_CATS; i++) {
#ifdef SOURCEQUAL_IGNORES_INDELS
          /* pretend it never happened: remove counts and ignore qualities */
          if (i==OP_INS || i==OP_DEL) {
               num_err_probs -= op_counts[i];
               continue;
          }
#endif
          if (i!=OP_MATCH) {
               num_non_matches += op_counts[i];
          }
          for (j=0; j<op_counts[i]; j++) {
               int qual;
               if (nonmatch_qual >= 0) {
                    qual = nonmatch_qual;
               } else {
                    qual = op_quals[i][j];
               }
               err_probs[err_prob_idx] = PHREDQUAL_TO_PROB(qual);
               /*LOG_FIXME("err_probs[%d] = %f (nonmatch_qual=%d op_quals[i=%d][j=%d]=%d)\n", err_prob_idx, err_probs[err_prob_idx], nonmatch_qual, i, j, op_quals[i][j]);*/
               err_prob_idx += 1;
          }
     }
     assert(err_prob_idx == num_err_probs);

     /*  need num_non_matches-1 */
     orig_num_non_matches = num_non_matches;
     if (num_non_matches>0) {
          num_non_matches -= 1;
     }
     if (0 == num_non_matches) {
#if 0
          src_qual = PROB_TO_PHREDQUAL_SAFE(0.0);
#else
          src_qual = PROB_TO_PHREDQUAL(LDBL_MIN);
#endif
          goto free_and_exit;
     }

#ifdef SOURCEQUAL_USES_PAIRS
     if ((b->core.flag&BAM_FPAIRED) && (b->core.flag&BAM_FPROPER_PAIR)) {
          double median_err = dbl_median(err_probs, num_err_probs);
          LOG_FIXME("%s\n", "SOURCEQUAL_USES_PAIRS: assuming perfect match for paire read");
          /* can't easily access info about read in pair. assume 
             perfect perfect match, using length and median prob of 
             current one */
          
          num_err_probs *= 2;
          if (NULL == (err_probs = realloc(err_probs, num_err_probs * sizeof(double)))) {
               fprintf(stderr, "FATAL: couldn't reallocate memory at %s:%s():%d\n",
                       __FILE__, __FUNCTION__, __LINE__);
               exit(1);
          }
          for (i=err_prob_idx-1; i<num_err_probs; i++) {
               err_probs[i] = median_err;
          }
          LOG_FIXME("median_err = %f\n", median_err);
     }
#endif

     /* src_prob: what's the prob of seeing n_mismatches-1 by chance,
      * given quals? or: how likely is this read from the genome.
      * 1-src_value = prob read is not from genome
      */

     /* sorting in theory should be numerically more stable and also
      * make poissbin faster */
     qsort(err_probs, num_err_probs, sizeof(double), dbl_cmp);
     probvec = poissbin(&unused_pval, err_probs,
                        num_err_probs, num_non_matches, 1.0, 0.05);
     /* need prob not pv */
     errno = 0;
     feclearexcept(FE_ALL_EXCEPT);
     src_prob = exp(probvec[num_non_matches-1]);
     if (errno || fetestexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW | FE_UNDERFLOW)) {
          if (src_prob < DBL_EPSILON) {
               src_prob = DBL_MIN;/* to zero but prevent actual 0 value */
          } else {
               src_prob = DBL_MAX; /* otherwise set to 1 which might pass filters */
          }
     }

     free(probvec);
     src_qual = PROB_TO_PHREDQUAL(1.0 - src_prob);


free_and_exit:
     for (i=0; i<NUM_OP_CATS; i++) {
          free(op_quals[i]);
     }
     free(op_quals);

     free(err_probs);

     /* if we wanted to use softening from precomputed stats then add
      * all non-matches up instead of using the matches */
#if 0
#define TRACE
#endif
#ifdef TRACE
     LOG_DEBUG("returning src_qual=%d (orig prob = %g) for cigar=%s num_err_probs=%d num_non_matches=%d(%d) @%d\n",
               src_qual, src_prob, cigar_str_from_bam(b), num_err_probs, num_non_matches, orig_num_non_matches, b->core.pos);
#endif
#undef TRACE

     return src_qual;
}
/* source_qual() */



/* not part of offical samtools/htslib API but part of samtools */
static int
mplp_func(void *data, bam1_t *b)
{
     mplp_aux_t *ma = (mplp_aux_t*)data;
     int ret, skip = 0;

     do {
          int has_ref;
          ret = ma->iter? sam_itr_next(ma->fp, ma->iter, b) : sam_read1(ma->fp, ma->h, b);
          if (ret < 0)
               break;

#ifdef TRACE
          LOG_DEBUG("Got read %s with flag %d\n", bam_get_qname(b), core.flag);
#endif
          if (b->core.tid < 0 || (b->core.flag&BAM_FUNMAP)) { /* exclude unmapped reads */
               skip = 1;
               continue;
          }
         
          if (b->core.flag & (BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP)) {
#ifdef TRACE
               LOG_DEBUG("%s BAM_DEF_MASK match\n", bam_get_qname(b));
#endif
               skip = 1; 
               continue;
          }
          if (ma->conf->bed) { /* test overlap */
               skip = !bed_overlap(ma->conf->bed, ma->h->target_name[b->core.tid], b->core.pos, bam_endpos(b));
               if (skip)
                    continue;
          }

          if (ma->conf->flag & MPLP_ILLUMINA13) {
               int i;
               uint8_t *qual = bam_get_qual(b);
               for (i = 0; i < b->core.l_qseq; ++i)
                    qual[i] = qual[i] > 31? qual[i] - 31 : 0;
          }
          has_ref = (ma->ref && ma->ref_id == b->core.tid)? 1 : 0;

          /* lofreq fix to original samtools routines which ensures that
           * the reads mapping to first position have a reference
           * attached as well and therefore baq, sq etc can be
           * applied */
          if (! has_ref && ma->conf->fai) {
               int ref_len = -1;
               ma->ref = faidx_fetch_seq(ma->conf->fai, ma->h->target_name[b->core.tid], 0, 0x7fffffff, &ref_len);
               if (!ma->ref) {
                    LOG_FATAL("Couldn't fetch sequence '%s'.\n", ma->h->target_name[b->core.tid]);
                    exit(1);/* FIXME just returning would just skip calls for this seq */

               } else {
                    strtoupper(ma->ref);/* safeguard */
                    ma->ref_id = b->core.tid;
                    has_ref = 1;
               }
          }

          skip = 0;
#if 0
          {
               fprintf(stderr, "before realn\n");
               samFile *fp = sam_open("-", "w");
               sam_hdr_write(fp, ma->h);
               sam_write1(fp, ma->h, b);
               fflush(stdout);
          }
#endif
          if (ma->conf->flag & MPLP_BAQ || ma->conf->flag & MPLP_IDAQ) {
               int baq_flag = ma->conf->flag & MPLP_BAQ ? 1 : 0;
               int baq_ext =  ma->conf->flag & MPLP_EXT_BAQ ? 1 : 0;
               int idaq_flag = ma->conf->flag & MPLP_IDAQ ? 1 : 0;

               if (! has_ref) {
                    LOG_FATAL("%s\n", "Can't compute BAQ or IDAQ without reference sequence");
                    exit(1);
               }
               if (baq_flag && ma->conf->flag & MPLP_REDO_BAQ) {
                    baq_flag = 2;
               }                    

               if (bam_prob_realn_core_ext(b, ma->ref, baq_flag, baq_ext, idaq_flag)) {
                    LOG_ERROR("bam_prob_realn_core() failed for %s\n", bam_get_qname(b));
               }

#if 0
               {
                    uint8_t *baq_aux = NULL;
                    baq_aux = bam_aux_get(b, BAQ_TAG);
                    if (! baq_aux) {
                         LOG_ERROR("bam_prob_realn_core() didn't report an error but %s is missing. Can happen on refskips etc\n", BAQ_TAG);
                    }

               }
#endif
          }

#if 0
          {
               fprintf(stdout, "after realn\n");
               samFile *fp = sam_open("-", "w");
               sam_hdr_write(fp, ma->h);
               sam_write1(fp, ma->h, b);
               fflush(stdout);
          }
      
#endif

        if (b->core.qual > ma->conf->max_mq) {
             b->core.qual = ma->conf->max_mq;
        } else if (b->core.qual < ma->conf->min_mq) {
             skip = 1;
        }
        /* never tried RELAXED_ORPHAN. most examples I saw where orphans matter evaluate to true under both conditions anyway */
#ifdef RELAXED_ORPHAN
        /* only orphan if mate is wrongly mapped (but not unmapped) */
        else if ((ma->conf->flag&MPLP_NO_ORPHAN) && (b->core.flag&BAM_FPAIRED) && !(b->core.flag&BAM_FMUNMAP)) {
#else
        /* orphan as in samtools: read is aligned but not as proper pair as defined by aligner */
        else if ((ma->conf->flag&MPLP_NO_ORPHAN) && (b->core.flag&BAM_FPAIRED) && !(b->core.flag&BAM_FPROPER_PAIR)) {
#endif
             skip = 1;
        }
    } while (skip);

    /* compute source qual if requested and have ref and attach as aux
     * to bam. only disadvantage of doing this here is that we don't
     * have BAQ info yet (only interesting if it's supposed to be used
     * instead of BQ) only have the ref but not the cons base.
     */
    if (ma->ref && ma->ref_id == b->core.tid && ma->conf->flag & MPLP_USE_SQ) {
         int sq = source_qual(b, ma->ref, ma->conf->def_nm_q,
                              ma->h->target_name[b->core.tid], DEFAULT_MIN_BQ/* FIXME could use->conf->min_bq which is set to a conservative 3 */);
         /* -1 indicates error or NA, but can't be stored as uint. hack is to use 0 instead */
         if (sq<0) {
              sq=0;
         }
         bam_aux_append(b, SRC_QUAL_TAG, 'i', sizeof(sq), (uint8_t*) &sq);
#if 0
         int sq2 = bam_aux2i(bam_aux_get(b, SRC_QUAL_TAG));
         LOG_WARN("sq=%d sq2=%d\n", sq, sq2);
#endif
    }

    return ret;
}

/* homopolymer run at (to the right of) current
                * position. if indels are not left aligned and current
                * position is already a homopolymer this will be taken
                * into account. mainly for filtering low af FP indel
                * at the beginning of poly-AT regions. A del GT>G
                * which is in the sequence context of GTTT will
                * receive an hrun value of 3. same for ins G>GT */
int
get_hrun(const int pos, const char *ref, const int ref_len)
{
     char c;
     int hrun=1;
     int i;

     /* to the right */
     i=pos+1;
     if (i>=ref_len) {
        return hrun;
     }   
     
     c=toupper(ref[i]);
     for (i=i+1; i<ref_len; i++) {
          /*LOG_DEBUG("to right: %c vs %c at %d\n", c, toupper(ref[i]), i+1);*/
          if (toupper(ref[i])==c) {
               hrun+=1;
          } else {
               break;
          }
     }

     /* extend to the left in case not left aligned */
     for (i=pos; i>=0; i--) {
          /*LOG_DEBUG("to left: %c vs %c at %d\n", c, toupper(ref[i]), i+1);*/
          if (toupper(ref[i])==c) {
               hrun+=1;
          } else {
               break;
          }
     }

     return hrun;
}

/* Press pileup info into one data-structure. plp_col members
 * allocated here. Called must free with plp_col_free();
 *
 * FIXME this used to be a convenience function and turned into a big
 * and slow monster. keeping copies of everything is inefficient and
 * in some cases an unnecessary waste of memory (e.g. SQ, BAQ or MQ if
 * not needed). Needs optimization.
 */
void compile_plp_col(plp_col_t *plp_col,
                 const bam_pileup1_t *plp, const int n_plp,
                 const mplp_conf_t *conf, const char *ref, const int pos,
                 const int ref_len, const char *target_name)
{
     int i;
     char ref_base;

     /* "base counts" minus error-probs before base-level filtering
      * for each base. temporary data-structure for cheaply determining
      * consensus which is saved in plp_col */
     double base_counts[NUM_NT4] = { 0 };
     /* sum of qualities for all non-indel events */
     int ins_nonevent_qual = 0, del_nonevent_qual = 0;

     /* computation of depth (after read-level *and* base-level filtering)
      * samtools-0.1.18/bam2depth.c:
      *   if (p->is_del || p->is_refskip) ++m;
      *   else if (bam_get_qual(p->b)[p->qpos] < bq) ++m
      * n_plp[i] - m
      */
     ref_base = (ref && pos < ref_len)? ref[pos] : 'N';

     plp_col_init(plp_col);
     plp_col->target = strdup(target_name);
     plp_col->pos = pos;
     plp_col->ref_base = ref_base;
     plp_col->coverage_plp = n_plp;  /* this is coverage as in the original mpileup,
                                    i.e. after read-level filtering */
     plp_col->num_bases = 0;
     plp_col->num_ign_indels = 0;
     plp_col->num_non_indels = 0;
     LOG_DEBUG("Processing %s:%d\n", plp_col->target, plp_col->pos+1);
     
     if (ref) {
          plp_col->hrun = get_hrun(pos, ref, ref_len);
     } else {
          plp_col->hrun = -1;
     }

     for (i = 0; i < n_plp; ++i) {
          /* inserted parts of pileup_seq() here.
           * logic there goes like this:
           *
           * if is_head: put ^
           *
           * if ! is_del: put c
           * else:        put * (or <> if refskip)
           *
           * if indel>0:   print + p[qpos +1...indel]
           * elif indel<0: print - p[qpos +1...indel]
           *
           * if is_tail: put $
           */
          const bam_pileup1_t *p = plp + i;
          int nt4;
          int mq=-1, bq, baq; /* phred scores */
          int iq = 0, dq = 0;
          int iaq = -1, daq = -1;
          int base_skip = 0; /* boolean */
          int sq = -1;
#ifdef USE_ALNERRPROF
          int aq = 0;
#endif
          uint8_t *bi = bam_aux_get(p->b, BI_TAG);
          uint8_t *bd = bam_aux_get(p->b, BD_TAG);
          uint8_t *ai = bam_aux_get(p->b, AI_TAG);
          uint8_t *ad = bam_aux_get(p->b, AD_TAG);
          uint8_t *baq_aux = NULL; /* full baq value (not offset as "BQ"!) */

#ifdef USE_OLD_AI_AD
          /* temporary fix preventing problems due to the fact that we changed AI AD to ai ad
           * to be deleted soon
           */
          if (! ai) {
               ai = bam_aux_get(p->b, "AI");
          }
          if (! ad) {
               ad = bam_aux_get(p->b, "AD");
          }
#endif
          if (conf->flag & MPLP_USE_SQ) {
               sq = bam_aux2i(bam_aux_get(p->b, SRC_QUAL_TAG)); /* lofreq internally computed on the fly */
          }

          if (conf->flag & MPLP_BAQ) {
               baq_aux = bam_aux_get(p->b, BAQ_TAG);
               /* should have been recomputed already */
               if (! baq_aux) {
                    if (! missing_baq_warning_printed) {
                         LOG_WARN("BAQ tag %s missing but BAQ was supposed to be used (at %s:%d; can happend if refskips are encountered)\n", BAQ_TAG, plp_col->target, plp_col->pos+1);
                         /*LOG_FATAL("%s\n", "Please pre-process your BAM file with lofreq alnqual first");*/
                         missing_baq_warning_printed = 1;
                    }
               } else {
                    baq_aux++; /* first char is type (same for bd and bi done below) */
               }
          }

#if 0
          LOG_FIXME("At %s:%d %c: p->is_del=%d p->is_refskip=%d p->indel=%d p->is_head=%d p->is_tail=%d\n",
                    plp_col->target,
                    plp_col->pos+1,
                    seq_nt16_str[bam_seqi(bam_get_seq(p->b), p->qpos)],
                    p->is_del, p->is_refskip, p->indel, p->is_head, p->is_tail);
#endif
          /* no need for check if mq is within user defined
           * limits. check was done in mplp_func */
          mq = p->b->core.qual;

          /* p->is_del mean there was a deletion (already printed before this column). */
          if (! p->is_del) {
               double count_incr;

               if (p->is_head) {
                    plp_col->num_heads += 1;
               }
               if (p->is_tail) {
                    plp_col->num_tails += 1;
               }

#if 0
               /* nt for printing */
               nt = seq_nt16_str[bam_seqi(bam_get_seq(p->b), p->qpos)];
               nt = bam_is_rev(p->b)? tolower(nt) : toupper(nt);
#endif

               /* nt4 for indexing */
               nt4 = seq_nt16_int[bam_seqi(bam_get_seq(p->b), p->qpos)];

               bq = bam_get_qual(p->b)[p->qpos];

               /* minimal base-call quality filtering. doing this here
                * will make all downstream analysis blind to filtering
                * and might skew AFs etc
                */
               if (bq < conf->min_plp_bq) {
                    base_skip = 1;
                    /* NOTE: all values used after check_indel need to be initialized before this goto */
                    goto check_indel; /* goto was easiest */
               }

               /* the following samtools' original code will correct
                * base-pairs down if they exceed the valid
                * sanger/phred limits. don't think it's wise to do this
                * automatically as this would indicate a problem with
                * the input and it's also unclear what the BAQ then means
                */
               if (bq > SANGER_PHRED_MAX) {
                    /* bq = SANGER_PHRED_MAX; /@ Sanger/Phred max */
                    LOG_WARN("Base quality above allowed maximum detected (%d > %d). Using max instead\n", bq, SANGER_PHRED_MAX, bam_get_qname(p->b));
                    bq = SANGER_PHRED_MAX;
               }
               PLP_COL_ADD_QUAL(& plp_col->base_quals[nt4], bq);

               if (baq_aux) {
                    baq = baq_aux[p->qpos]-33;
                    PLP_COL_ADD_QUAL(& plp_col->baq_quals[nt4], baq);
               } else if (conf->flag & MPLP_BAQ)  {
                    /* baq was enabled but failed. set to -1 */
                    PLP_COL_ADD_QUAL(& plp_col->baq_quals[nt4], -1);
               }

               /* samtools check to detect Sanger max value: problem
                * is that an MQ Phred of 255 means NA according to the
                * samtools spec (needed below). This however is not
                * detectable if the following original samtools line
                * gets executed, which is why we remove it:
                * if (mq > 126) mq = 126;
                */
               PLP_COL_ADD_QUAL(& plp_col->map_quals[nt4], mq);

               if (conf->flag & MPLP_USE_SQ) {
                    PLP_COL_ADD_QUAL(& plp_col->source_quals[nt4], sq);
               }
#ifdef USE_ALNERRPROF
               if (alnerrprof) {
                    int tid = p->b->core.tid;
                    assert(tid < alnerrprof->num_targets);
                    if (alnerrprof->prop_len[tid] > p->qpos) {
                         aq = PROB_TO_PHREDQUAL_SAFE(alnerrprof->props[tid][p->qpos]);
                         PLP_COL_ADD_QUAL(& plp_col->alnerr_qual[nt4], aq);
                    } else {
                         LOG_ERROR("alnerror for tid=%d too small for qpos=%d. Setting to 0\n", tid, p->qpos+1);
                         PLP_COL_ADD_QUAL(& plp_col->alnerr_qual[nt4], PROB_TO_PHREDQUAL(LDBL_MIN));
                    }
               }
               /* don't add anything. keep empty */
#endif


#ifdef MERGEQ_FOR_CONS_CALL
#ifdef USE_ALNERRPROF
               count_incr = 1.0 - merge_srcq_baseq_mapq_and_alnq(sq, bq, mq, aq);
#else
               count_incr = 1.0 - merge_srcq_baseq_and_mapq(sq, bq, mq);
#endif
#else
               count_incr = 1.0 - PHREDQUAL_TO_PROB(bq);
#endif

               /* FIXME this can't be the proper way to handle cases where count_incr = 0.0 because one of the values is 0? */
               if (count_incr == 0.0) {
                    count_incr = DBL_MIN;
               }

               base_counts[nt4] += count_incr;
               if (bam_is_rev(p->b)) {
                    plp_col->rv_counts[nt4] += 1;
               } else {
                    plp_col->fw_counts[nt4] += 1;
               }

          } /* ! p->is_del */
          /* else {deletion (already printed before this column), i.e. we got physical coverage (if no terminal indels are allowed} */

check_indel:

          /* for post read- and base-level coverage. FIXME review */
          if (! (p->is_del || p->is_refskip || 1 == base_skip)) {/* FIXME also use p->indel? */
               plp_col->num_bases += 1;
          }

          if (bi) {
               char *t = (char*)(bi+1); /* 1 is type */
#if 0
               int j;
               printf("At %d qpos %d: %s=%s", plp_col->pos+1, p->qpos+1, BI_TAG, t);
               for (j = 0; j < p->indel; ++j) {
                    printf(" %c:%d-%d-%d", t[p->qpos+j], t[p->qpos+j-1]-33, t[p->qpos+j]-33, t[p->qpos+j+1]-33);
               }
               printf("\n");
#endif
               /* adding 1 value representing whole insert */
               iq = t[p->qpos] - 33;
          } /* else default to 0 */

          if (bd) {
               char *t = (char*)(bd+1);  /* 1 is type */
#if 0
               int j;
               printf("At %d qpos %d: %s=%sx", plp_col->pos+1, p->qpos+1, BD_TAG, t);
               for (j = 0; j < p->indel; ++j) {
                    printf(" %c:%d-%d-%d", t[p->qpos+j], t[p->qpos+j-1]-33, t[p->qpos+j]-33, t[p->qpos+j+1]-33);
               }
               printf("\n");
#endif
               /* adding 1 value representing whole del */
               dq = t[p->qpos] - 33;
#ifdef PACBIO_REALN_HRUNDQ7
               /* FIXME temp artifically decreasing pacbio del quals in hruns */
               if (plp_col->hrun>1) {
                    if (dq>7) dq=7;
               }
#elif PACBIO_REALN
               /* FIXME temp artifically decreasing pacbio del quals in hruns */
               if (dq>=10) dq-=10;
#endif
          } /* else default to 0 */


          if (iq < conf->min_plp_idq || dq < conf->min_plp_idq) {
               /* LOG_DEBUG("iq=%d < conf->min_plp_idq=%d || dq=%d < conf->min_plp_idq=%d\n", iq, conf->min_plp_idq, dq, conf->min_plp_idq); */
               if (p->indel != 0 || p->is_del != 0) {
                  plp_col->num_ign_indels += 1;
               }
          } else {

               if (p->indel != 0) {
                    /* insertion (+)
                     */
                    if (p->indel > 0) {
                         char *ins_seq;
                         int j;

                         if (ai) {
                              char *a = (char*)(ai+1);
                              iaq = a[p->qpos] - 33;
                              plp_col->has_indel_aqs = 1;
                         }

                         plp_col->num_ins += 1;
                         plp_col->sum_ins += p->indel;

                         if ((ins_seq = malloc((p->indel+1) * sizeof(char)))==NULL) {
                              LOG_FATAL("%s\n", "Memory allocation failed");
                              exit(1);
                         }

                         /* get inserted sequence */
                         for (j = 1; j <= p->indel; ++j) {
                              int c = seq_nt16_str[bam_seqi(bam_get_seq(p->b), p->qpos+j)];
                              ins_seq[j-1] = toupper(c);
                         }
                         ins_seq[j-1] = '\0';


                         /*LOG_DEBUG("Insertion of %s at %d with iq %d iaq %d\n", ins_seq, pos, iq, iaq);*/
                         add_ins_sequence(&plp_col->ins_event_counts,
                              ins_seq, iq, iaq, mq, sq,
                              bam_is_rev(p->b)? 1: 0);

                         PLP_COL_ADD_QUAL(& plp_col->del_quals, dq);
                         PLP_COL_ADD_QUAL(& plp_col->del_map_quals, mq);
                         PLP_COL_ADD_QUAL(& plp_col->del_source_quals, sq);
                         del_nonevent_qual += dq;
                         if (bam_is_rev(p->b)) {
                              plp_col->non_del_fw_rv[1] += 1;
                         } else {
                              plp_col->non_del_fw_rv[0] += 1;
                         }
                         free(ins_seq);

                    /* deletion (-)
                     */
                    } else if (p->indel < 0) {
                         /* get deleted sequence */
                         char *del_seq;
                         int j;

                         if (ad) {
                              char *a = (char*)(ad+1);
                              daq = a[p->qpos] - 33;
                              plp_col->has_indel_aqs = 1;
                         }

                         plp_col->num_dels += 1;
                         plp_col->sum_dels -= p->indel;

                         if ((del_seq = malloc(((-p->indel)+1) * sizeof(char)))==NULL) {
                              LOG_FATAL("%s\n", "Memory allocation failed");
                              exit(1);
                         }

                         for (j = 1; j <= -p->indel; ++j) {
                              int c =  (ref && (int)pos+j < ref_len)? ref[pos+j] : 'N';
                              del_seq[j-1] = toupper(c);
                         }
                         del_seq[j-1] = '\0';
                         /*LOG_DEBUG("Deletion of %s at %d with dq %d daq %d\n", del_seq, pos, dq, daq);*/
#ifdef PACBIO_SUPPRESS_1BASE_DEL
                         /*LOG_FIXME("Deletion of %s at %d with dq %d daq %d\n", del_seq, pos, dq, daq);*/
                         if (strlen(del_seq)==1) {
                              if (del_seq[0]=='G' || del_seq[0]=='C') {
                                   dq -= 10;
                              }
                              if (del_seq[0]=='A' || del_seq[0]=='T') {
                                   dq -= 5;
                              }
                              if (dq<0) {
                                   dq=0;
                              }
                         }
#endif
                         add_del_sequence(&plp_col->del_event_counts,
                              del_seq, dq, daq, mq, sq,
                              bam_is_rev(p->b)? 1: 0);
                         PLP_COL_ADD_QUAL(& plp_col->ins_quals, iq);
                         PLP_COL_ADD_QUAL(& plp_col->ins_map_quals, mq);
                         PLP_COL_ADD_QUAL(& plp_col->ins_source_quals, sq);
                         ins_nonevent_qual += iq;
                         if (bam_is_rev(p->b)) {
                              plp_col->non_ins_fw_rv[1] += 1;
                         } else {
                              plp_col->non_ins_fw_rv[0] += 1;
                         }
                         free(del_seq);
                    }

               } else { /* if (p->indel != 0) ... */
                    plp_col->num_non_indels += 1;
                    /* neither deletion, nor insertion. need the qualities anyway */
                    PLP_COL_ADD_QUAL(& plp_col->ins_quals, iq);
                    PLP_COL_ADD_QUAL(& plp_col->ins_map_quals, mq);
                    ins_nonevent_qual += iq;
                    if (bam_is_rev(p->b)) {
                         plp_col->non_ins_fw_rv[1] += 1;
                    } else {
                         plp_col->non_ins_fw_rv[0] += 1;
                    }

                    /*LOG_DEBUG("Neither deletion nor insertion. Adding iq=%d dq=%d\n", iq, dq);*/
                    PLP_COL_ADD_QUAL(& plp_col->del_quals, dq);
                    PLP_COL_ADD_QUAL(& plp_col->del_map_quals, mq);
                    del_nonevent_qual += dq;
                    if (bam_is_rev(p->b)) {
                         plp_col->non_del_fw_rv[1] += 1;
                    } else {
                         plp_col->non_del_fw_rv[0] += 1;
                    }
               }
          }

     }  /* end: for (i = 0; i < n_plp; ++i) { */


     /* ****************** FINDING CONSENSUS **************** */
     /* consensus is saved as a char array starting with '+' or '-'
      * if the consensus is an insertion or deletion. there is an
      * insertion event if the sum of qualities for that insertion event
      * is greater than the sum of qualities for all non-insertion events.
      * there is a deletion event if the sum of qualities for that
      * deletion event is greater than the sum of qualities for all non-deletion
      * events. otherwise, the consensus base is not an indel and is given
      * by the nucleotide with the greatest sum of qualities.
      *
      *
      * YHT 2/10/14: """the idea is to find the event which has the highest probability of
      * occurring the number of times it was observed. For example, if I see
      * 2 +A events at a position, the probability that those 2 events are real
      * and occurred together is (1 - the probability that the +A event occurred
      * due to error for the first supporting read)*(1 - the probability that the
      * +A event occurred due to error for the second supporting read). I keep
      * track of this for each insertion event, for e.g. +AA, +T etc. at that
      * position. I guess in theory this should incorporate errors from other sources.
      *
      * I also find the probability of seeing n number of non-insertions at that
      * position, which is the product of (1 - the probability of seeing an insertion
      * error at that position for each of those reads).
      *
      * I then compare the probabilities of each of these events. It turns out that
      * comparing the products of (1 - error probability) is the same as comparing
      * the log sum of the qualities, because qualities are the negative log of the
      * error probabilities."""
      *
      * FIXME: check consensus indel against minimum consensus quality
      * FIXME: merge indel qualities when determining consensus indel event
      * FIXME(AW): why are we using max qualities here and not errprob corrected counts?
      */

     ins_event *ins_it, *ins_it_tmp;
     char *ins_maxevent_key = NULL;
     int ins_maxevent_qual = 0;
     HASH_ITER(hh_ins, plp_col->ins_event_counts, ins_it, ins_it_tmp) {
          if (ins_it->cons_quals > ins_maxevent_qual) {
               ins_maxevent_key = ins_it->key;
               ins_maxevent_qual = ins_it->cons_quals;
          }
     }
     del_event *del_it, *del_it_tmp;
     char *del_maxevent_key = NULL;
     int del_maxevent_qual = 0;
     HASH_ITER(hh_del, plp_col->del_event_counts, del_it, del_it_tmp) {
          if (del_it->cons_quals > del_maxevent_qual) {
               del_maxevent_key = del_it->key;
               del_maxevent_qual = del_it->cons_quals;
          }
     }

     /* LOG_DEBUG("ins_maxevent_qual:%d ins_nonevent_qual:%d "
               "del_maxevent_qual:%d del_nonevent_qual:%d\n",
               ins_maxevent_qual, ins_nonevent_qual,
               del_maxevent_qual, del_nonevent_qual); */

     if (!(ins_maxevent_qual > ins_nonevent_qual) &&
         !(del_maxevent_qual > del_nonevent_qual)) {
          /* determine consensus from 'counts'. will never produce N on tie  */
          plp_col->cons_base[0] = bam_nt4_rev_table[
               argmax_d(base_counts, NUM_NT4)];
          plp_col->cons_base[1] = '\0';
     } else if (ins_maxevent_qual > ins_nonevent_qual) {  // consensus insertion
          /* LOG_DEBUG("cons ins: ins_maxevent_qual=%d > ins_nonevent_qual=%d\n", ins_maxevent_qual, ins_nonevent_qual); */
          plp_col->cons_base[0] = '+';
          strcpy(plp_col->cons_base+1, ins_maxevent_key);
     } else if (del_maxevent_qual > del_nonevent_qual) { // consensus deletion
          /* LOG_DEBUG("cons del: del_maxevent_qual=%d > del_nonevent_qual=%d\n", del_maxevent_qual, del_nonevent_qual); */
          plp_col->cons_base[0] = '-';
          strcpy(plp_col->cons_base+1, del_maxevent_key);
     } else {
          LOG_FATAL("internal error...");
          exit(1);
     }

     if (debug) {
          plp_col_debug_print(plp_col, stderr);
     }
#if 0
     plp_col_mpileup_print(plp_col, stdout);
#endif

     for (i = 0; i < NUM_NT4; ++i) {
          assert(plp_col->fw_counts[i] + plp_col->rv_counts[i] == plp_col->base_quals[i].n);
          assert(plp_col->base_quals[i].n == plp_col->baq_quals[i].n);
          assert(plp_col->base_quals[i].n == plp_col->map_quals[i].n);
          assert(plp_col->map_quals[i].n == plp_col->source_quals[i].n);
     }
}
/* compile_plp_col() */



/* not part of offical samtools/htslib API but part of samtools */
int
mpileup(const mplp_conf_t *mplp_conf,
        void (*plp_proc_func)(const plp_col_t*, void*),
        void *plp_proc_conf,
        const int n, const char **fn)
{
    mplp_aux_t **data;
    int i, tid, pos, *n_plp, tid0 = -1, beg0 = 0, end0 = 1u<<29, ref_len = -1, ref_tid = -1, max_depth;
    const bam_pileup1_t **plp;
    bam_mplp_t iter;
    bam_hdr_t *h = 0;
    char *ref;
    kstring_t buf;
    long long int plp_counter = 0; /* note: some cols are simply skipped */

    /* paranoid exit. n only allowed to be one in our case (not much
     * of an *m*pileup, I know...) */
    if (1 != n) {
         fprintf(stderr, "FATAL(%s:%s): need exactly one BAM files as input (got %d)\n",
                 __FILE__, __FUNCTION__, n);
         for (i=0; i<n; i++) {
              fprintf(stderr, "%s\n", fn[i]);
         }
         return 1;
    }

    memset(&buf, 0, sizeof(kstring_t));
    data = calloc(n, sizeof(mplp_aux_t*));
    plp = calloc(n, sizeof(bam_pileup1_t*));
    n_plp = calloc(n, sizeof(int));


    /* read the header and initialize data
     *
     * note: most of this is overkill since it deals with multiple bam
     * files, whereas we allow only one. however, if we keep it close
     * to the original source then a diff against future versions of
     * samtools is easier
     *
     */
    for (i = 0; i < n; ++i) {
        bam_hdr_t *h_tmp;
        if (0 != strcmp(fn[i], "-")) {
          if (! file_exists(fn[i])) {
            fprintf(stderr, "File '%s' does not exist. Exiting...\n", fn[i]);
            exit(1);
          }
        }
        data[i] = calloc(1, sizeof(mplp_aux_t));
        data[i]->fp = sam_open(fn[i], "r");
        data[i]->conf = mplp_conf;
        h_tmp = sam_hdr_read(data[i]->fp);
        if ( !h_tmp ) {
             fprintf(stderr,"[%s] fail to read the header of %s\n", __func__, fn[i]);
             exit(1);
        }
        data[i]->h = i? h : h_tmp; /* for i==0, "h" has not been set yet */

        if (mplp_conf->reg) {
            hts_idx_t *idx;
            idx = sam_index_load(data[i]->fp, fn[i]);
            if (idx == 0) {
                fprintf(stderr, "[%s] fail to load index for %d-th input.\n", __func__, i+1);
                exit(1);
            }
            if ((data[i]->iter = sam_itr_querys(idx, h_tmp, mplp_conf->reg)) == NULL) {
                fprintf(stderr, "[%s] malformatted region or wrong seqname for %d-th input.\n", __func__, i+1);
                exit(1);
            }
            if (i == 0) tid0 = data[i]->iter->tid, beg0 = data[i]->iter->beg, end0 = data[i]->iter->end;
            hts_idx_destroy(idx);
        }
        if (i == 0) {
             h = h_tmp;
        } else {
            bam_hdr_destroy(h_tmp);
        }
    }
    LOG_DEBUG("%s\n", "BAM header initialized");
    if (debug) {
         for (i=0; i < h->n_targets; i++) {
              LOG_DEBUG("BAM header target #%d: name=%s len=%d\n", i, h->target_name[i], h->target_len[i]);
         }
    }
    if (tid0 >= 0 && mplp_conf->fai) { /* region is set */
         ref = faidx_fetch_seq(mplp_conf->fai, h->target_name[tid0], 0, 0x7fffffff, &ref_len);
         if (NULL == ref || h->target_len[tid0] != ref_len) {
              LOG_FATAL("Reference fasta file doesn't seem to contain the right sequence(s) for this BAM file. (mismatch for seq %s listed in BAM header)\n", h->target_name[tid0]);
              return -1;
         }
         strtoupper(ref);/* safeguard */
         ref_tid = tid0;
         for (i = 0; i < n; ++i) data[i]->ref = ref, data[i]->ref_id = tid0;
    } else {
         ref_tid = -1;
         ref = 0;
    }
    iter = bam_mplp_init(n, mplp_func, (void**)data);
    max_depth = mplp_conf->max_depth;
    bam_mplp_set_maxcnt(iter, max_depth);

#ifdef USE_ALNERRPROF
    if (mplp_conf->alnerrprof_file) {
         alnerrprof = calloc(1, sizeof(alnerrprof_t));
         if (parse_alnerrprof_statsfile(alnerrprof, mplp_conf->alnerrprof_file, h)) {
              LOG_FATAL("parse_errprof_statsfile() on %s failed\n", mplp_conf->alnerrprof_file);
              exit(1);
         }
         normalize_alnerrprof(alnerrprof);
    }
#endif

    LOG_DEBUG("%s\n", "Starting pileup loop");
    while (bam_mplp_auto(iter, &tid, &pos, n_plp, plp) > 0) {
        plp_col_t plp_col;
        int i=0; /* NOTE: mpileup originally iterated over n */

        if (mplp_conf->reg && (pos < beg0 || pos >= end0))
             continue; /* out of the region requested */
        if (mplp_conf->bed && tid >= 0 && !bed_overlap(mplp_conf->bed, h->target_name[tid], pos, pos+1))
             continue;
        if (tid != ref_tid) {
            free(ref); ref = 0;
            if (mplp_conf->fai) {
                 ref = faidx_fetch_seq(mplp_conf->fai, h->target_name[tid], 0, 0x7fffffff, &ref_len);
                 if (NULL == ref || h->target_len[tid] != ref_len) {
                      LOG_DEBUG("ref %s at %p h->target_len[tid]=%d ref_len=%d\n", h->target_name[tid], ref, h->target_name[tid], ref_len)
                      LOG_FATAL("Reference fasta file doesn't seem to contain the right sequence(s) for this BAM file. (mismatch for seq %s listed in BAM header).\n", h->target_name[tid]);
                      return -1;
                 }
                 strtoupper(ref);/* safeguard */
                 LOG_DEBUG("%s\n", "sequence fetched");
            }
            for (i = 0; i < n; ++i)  {
                 data[i]->ref = ref, data[i]->ref_id = tid;
            }
            ref_tid = tid;
        }
        i=0; /* i is 1 for first pos which is a bug due to the removal
              * of one of the loops, so reset here */

        plp_counter += 1;
        if (1 == plp_counter%100000) {
             LOG_VERBOSE("Alive and happily crunching away on pos"
                         " %d of %s...\n", pos+1, h->target_name[tid]);
        }

        compile_plp_col(&plp_col, plp[i], n_plp[i], mplp_conf,
                        ref, pos, ref_len, h->target_name[tid]);

        (*plp_proc_func)(& plp_col, plp_proc_conf);

        plp_col_free(& plp_col);

    } /* while bam_mplp_auto */

#ifdef USE_ALNERRPROF
    if (alnerrprof) {
         free_alnerrprof(alnerrprof);
         free(alnerrprof);
    }
#endif
    free(buf.s);
    bam_mplp_destroy(iter);
    bam_hdr_destroy(h);
    for (i = 0; i < n; ++i) {
        sam_close(data[i]->fp);
        if (data[i]->iter) bam_itr_destroy(data[i]->iter);
        free(data[i]);
    }
    free(data); free(plp); free(ref); free(n_plp);
    return 0;
}
/* mpileup() */
