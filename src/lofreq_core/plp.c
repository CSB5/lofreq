/* -*- c-file-style: "k&r"; indent-tabs-mode: nil; -*-
 *
 * This file is partially based on samtools' bam_plcmd.c and very
 * likely needs an update whenever samtools/libbam is updated
 *
 * FIXME missing license
 *
 */
#include <ctype.h>
#include <assert.h>

#include "faidx.h"
#include "kstring.h"
#include "sam.h"
#include "log.h"
#include "plp.h"
#include "vcf.h"
#include "samutils.h"
#include "snpcaller.h"


/* from bedidx.c */
void *bed_read(const char *fn);
void bed_destroy(void *_h);
int bed_overlap(const void *_h, const char *chr, int beg, int end);

/* Using XS, XG and ZS already used by others. So use ZG for
 * source qual. Not modyfying BAM anyway, so it wouldn't
 * matter if we overwrite an existing tag. */
#define SRC_QUAL_TAG "ZG"

const char *bam_nt4_rev_table = "ACGTN";

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


#ifdef USE_ALNERRPROF
static alnerrprof_t *alnerrprof = NULL;
#endif

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

    const int grow_by_size = 1000;

    p->target =  NULL;
    p->pos = -INT_MAX;
    p->ref_base = '\0';
    p->cons_base = 'N';
    p->coverage = -INT_MAX;
    for (i=0; i<NUM_NT4; i++) {
         int_varray_init(& p->base_quals[i], grow_by_size);
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
    int_varray_init(& p->ins_quals, 0);
    p->num_dels = p->sum_dels = 0;
    int_varray_init(& p->del_quals, 0);
}


void
plp_col_free(plp_col_t *p) {
    int i;

    free(p->target);
    for (i=0; i<NUM_NT4; i++) {
         int_varray_free(& p->base_quals[i]);
         int_varray_free(& p->map_quals[i]);
         int_varray_free(& p->source_quals[i]);
#ifdef USE_ALNERRPROF
         int_varray_free(& p->alnerr_qual[i]);
#endif
    }
    int_varray_free(& p->ins_quals);
    int_varray_free(& p->del_quals);
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
void
plp_col_mpileup_print(const plp_col_t *p, FILE *stream)
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
/* plp_col_mpileup_print() */




void
dump_mplp_conf(const mplp_conf_t *c, FILE *stream) 
{
     fprintf(stream, "mplp options\n");
     fprintf(stream, "  max_mq       = %d\n", c->max_mq);
     fprintf(stream, "  min_mq       = %d\n", c->min_mq);
     fprintf(stream, "  flag         = %d\n", c->flag);

     fprintf(stream, "  flag & MPLP_NO_ORPHAN  = %d\n", c->flag&MPLP_NO_ORPHAN ? 1:0);
     fprintf(stream, "  flag & MPLP_REALN      = %d\n", c->flag&MPLP_REALN ? 1:0);
     fprintf(stream, "  flag & MPLP_USE_SQ     = %d\n", c->flag&MPLP_USE_SQ ? 1:0);
     fprintf(stream, "  flag & MPLP_REDO_BAQ    = %d\n", c->flag&MPLP_REDO_BAQ ? 1:0);
     fprintf(stream, "  flag & MPLP_ILLUMINA13 = %d\n", c->flag&MPLP_ILLUMINA13 ? 1:0);
     
     fprintf(stream, "  capQ_thres   = %d\n", c->capQ_thres);
     fprintf(stream, "  max_depth    = %d\n", c->max_depth);
     fprintf(stream, "  min_bq       = %d\n", c->min_bq);
     fprintf(stream, "  def_nm_q     = %d\n", c->def_nm_q);
     fprintf(stream, "  reg          = %s\n", c->reg);
     fprintf(stream, "  fa           = %p\n", c->fa);
     /*fprintf(stream, "  fai          = %p\n", c->fai);*/
     fprintf(stream, "  bed          = %p\n", c->bed);
     fprintf(stream, "  cmdline      = %s\n", c->cmdline);
}
/* dump_mplp_conf() */



/* FIXME get rid of function in future */
static inline int
printw(int c, FILE *fp)
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


#ifdef USE_SOURCEQUAL

static var_hash_t *source_qual_ign_vars_hash = NULL; /* must be declared NULL ! */


int 
var_in_ign_list(var_t *var) {
     char *key = NULL;
     var_hash_t *match = NULL;

     /* using key_simple i.e. chrom and pos only */
     vcf_var_key_simple(&key, var);
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
         if (-1 == rc) {
              LOG_FATAL("%s\n", "Parsing error while parsing 2nd vcf-file");
              exit(1);
         }
         if (1 == rc) {/* EOF */
              free(var);
              break;
         }
         num_total_vars += 1;

         if (! read_only_passed || VCF_VAR_PASSES(var)) {
              var_hash_t *match = NULL;

              if (bed && ! bed_overlap(bed, var->chrom, var->pos, var->pos+1)) {
                   continue;
              }
              /* using key_simple i.e. chrom and pos only */
              vcf_var_key_simple(&key, var);

              HASH_FIND_STR(source_qual_ign_vars_hash, key, match);
              if (match) {
                   LOG_DEBUG("Already got a variant match for key '%s'. Will keep the old one.\n", key);
                   free(var);
                   continue;
              }

              var_hash_add(& source_qual_ign_vars_hash, key, var);
         }

#ifdef TRACE
         LOG_DEBUG("Adding %s\n", key);
#endif
    }

    num_kept_vars = HASH_COUNT(source_qual_ign_vars_hash);
    if (num_kept_vars) {
         LOG_VERBOSE("Kept %d variants (of a total of %d) to ignore from %s\n",
                     num_kept_vars, num_total_vars, vcf_path);
    } else {
         LOG_WARN("None of the %d variants in %s were kept\n", 
                  num_total_vars, vcf_path);
    }
    vcf_file_close(& vcf_file);

    return 0;
}


#endif


/* Estimate as to how likely it is that this read, given the mapping,
 * comes from this reference genome. P(r not from g|mapping) = 1 - P(r
 * from g). 
 * 
 * The overall idea was to use something as follows:
 * PJ = PM  +  (1-PM) * PS  +  (1-PM) * (1-PS) * PB, where
 * PJ = joined error probability
 * PM = mapping err.prob.
 * PS = source/genome err.prob.
 * PB = base err.prob.
 * 
 * In theory PS should go first but the rest is hard to compute then.
 * Using PM things get tractable and it intrinsically takes care of
 * PS.
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
 * Returns -1 on error. otherwise prob to see the observed number of
 * mismatches-1
 *
 */
int
source_qual(const bam1_t *b, const char *ref, const int nonmatch_qual, char *target)
{
     int op_counts[NUM_OP_CATS];
     int **op_quals = NULL;

     double *probvec = NULL;
     int num_non_matches; /* incl indels */
     int orig_num_non_matches = -1;
     double *err_probs = NULL; /* error probs (qualities) passed down to snpcaller */
     int num_err_probs; /* #elements in err_probs */

     double unused_pval;
     int src_qual = 255;
     double src_prob = -1; /* prob of this read coming from genome */
     int err_prob_idx;
     int i, j;

     /* alloc op_quals
      */
     if (NULL == (op_quals = malloc(NUM_OP_CATS * sizeof(int *)))) {
          fprintf(stderr, "FATAL: couldn't allocate memory at %s:%s():%d\n",
                  __FILE__, __FUNCTION__, __LINE__);
          exit(1);
     }
     for (i=0; i<NUM_OP_CATS; i++) {
          if (NULL == (op_quals[i] = malloc(MAX_READ_LEN * sizeof(int)))) {
               fprintf(stderr, "FATAL: couldn't allocate memory at %s:%s():%d\n",
                       __FILE__, __FUNCTION__, __LINE__);
               free(op_quals);
               exit(1);
          }
     }
          
     /* count match operations and get qualities for them
      */
/*     LOG_FIXME("%s\n", "Don't know ref name in count_cigar_ops which would be needed as hash key");*/
     num_err_probs = count_cigar_ops(op_counts, op_quals, b, ref, -1, target);
     if (-1 == num_err_probs) {
          LOG_WARN("%s\n", "count_cigar_ops failed on read"); /* FIXME print read */
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
         if (i != OP_MATCH) {
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
          src_qual = PROB_TO_PHREDQUAL(0.0);
          goto free_and_exit;
     }

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
     src_prob = exp(probvec[num_non_matches-1]);
     free(probvec);
     src_qual = PROB_TO_PHREDQUAL(1.0 - src_prob);


free_and_exit:
     for (i=0; i<NUM_OP_CATS; i++) {
          free(op_quals[i]);
     } 
     free(op_quals);

     free(err_probs);

     /* if we wanted to use softening from precomputed stats then add all non-matches up instead of using the matches */
#ifdef TRACE
     LOG_DEBUG("returning src_qual=%d (orig prob = %g) for cigar=%s num_err_probs=%d num_non_matches=%d(%d) @%d\n", 
               src_qual, src_prob, cigar_str_from_bam(b), num_err_probs, num_non_matches, orig_num_non_matches, b->core.pos);
#endif
     return src_qual;
}
/* source_qual() */



/* modelled after samtools FIXME.c */
static int
mplp_func(void *data, bam1_t *b)
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
        if (has_ref && (ma->conf->flag & MPLP_REALN)) {
             bam_prob_realn_core(b, ma->ref, (ma->conf->flag & MPLP_REDO_BAQ)? 7 : 3);
        }
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

#ifdef USE_SOURCEQUAL
    /* compute source qual if requested and have ref and attach as aux to bam */
    if (ma->ref && ma->ref_id == b->core.tid && ma->conf->flag & MPLP_USE_SQ) {
         int sq = source_qual(b, ma->ref, ma->conf->def_nm_q, ma->h->target_name[b->core.tid]);
          /* see bam_md.c for examples of bam_aux_append()
          * FIXME only allows us to store values as uint8_t i.e. a byte, i.e. 255 is max (that's also why len 4)
          */
         if (sq>255) {
              sq=255;
         }
         bam_aux_append(b, SRC_QUAL_TAG, 'i', 4, (uint8_t*) &sq);     
    }
#endif
    return ret;
}



/* Convenience function to press pileup info into one easy to handle
 * data-structure. plp_col members allocated here. Called must free
 * with plp_col_free(plp_col);
 *
 * FIXME this used to be a convenience function and turned into a big
 * and slow monster. keeping copies of everything is inefficient.
 * const bam_pileup1_t *plp is (almost) all we need for snv calling,
 * so get rid of this function in the future (which will however break
 * the subcommand plpsummmary)
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

     /* computation of depth (after read-level *and* base-level filtering)
      * samtools-0.1.18/bam2depth.c: 
      *   if (p->is_del || p->is_refskip) ++m; 
      *   else if (bam1_qual(p->b)[p->qpos] < bq) ++m
      * n_plp[i] - m
      */
     int num_skips = 0;

     ref_base = (ref && pos < ref_len)? ref[pos] : 'N';
     
     plp_col_init(plp_col);
     plp_col->target = strdup(target_name);
     plp_col->pos = pos;
     plp_col->ref_base = toupper(ref_base);
     plp_col->coverage = n_plp;  /* this is coverage as in the original mpileup, 
                                   i.e. after read-level filtering */
     LOG_DEBUG("Processing %s:%d\n", plp_col->target, plp_col->pos+1);
     
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
#ifdef USE_SOURCEQUAL
          int sq = -1;
          if (conf->flag & MPLP_USE_SQ) {
               sq = bam_aux2i(bam_aux_get(p->b, SRC_QUAL_TAG)); /* lofreq internally computed on the fly */
          }
#endif

#if 0
          LOG_FIXME("At %s:%d %c: p->is_del=%d p->is_refskip=%d p->indel=%d p->is_head=%d p->is_tail=%d\n", 
                    plp_col->target, 
                    plp_col->pos+1,
                    bam_nt16_rev_table[bam1_seqi(bam1_seq(p->b), p->qpos)],
                    p->is_del, p->is_refskip, p->indel, p->is_head, p->is_tail);
#endif

          if (! p->is_del) {
               if (p->is_head) {
                    plp_col->num_heads += 1;
               }
               if (p->is_tail) {
                    plp_col->num_tails += 1;
               }

#if 0
               /* nt for printing */
               nt = bam_nt16_rev_table[bam1_seqi(bam1_seq(p->b), p->qpos)];
               nt = bam1_strand(p->b)? tolower(nt) : toupper(nt);
#endif

               /* nt4 for indexing */
               nt4 = bam_nt16_nt4_table[bam1_seqi(bam1_seq(p->b), p->qpos)];                           

               bq = bam1_qual(p->b)[p->qpos];

               /* minimal base-call quality filtering
                */
               if (bq < conf->min_bq) {
                    base_skip = 1;
                    goto check_indel; /* goto was easiest */
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

               PLP_COL_ADD_QUAL(& plp_col->base_quals[nt4], bq);
               PLP_COL_ADD_QUAL(& plp_col->map_quals[nt4], mq);
#ifdef USE_SOURCEQUAL
               if (conf->flag & MPLP_USE_SQ) {
                    PLP_COL_ADD_QUAL(& plp_col->source_quals[nt4], sq);
               }
#endif
#ifdef USE_ALNERRPROF
               if (alnerrprof) {
                    int tid = p->b->core.tid;
                    assert(tid < alnerrprof->num_targets);
                    if (alnerrprof->prop_len[tid] > p->qpos) {
                         int q;
                         q = PROB_TO_PHREDQUAL(alnerrprof->props[tid][p->qpos]);
                         PLP_COL_ADD_QUAL(& plp_col->alnerr_qual[nt4], q);
                    } else {
                         LOG_ERROR("alnerror for tid=%d too small for qpos=%d. Setting to 0\n", tid, p->qpos+1);
                         PLP_COL_ADD_QUAL(& plp_col->alnerr_qual[nt4], PROB_TO_PHREDQUAL(0.0));
                    }
               }
               /* don't add anything. keep empty */
#endif
               if (bam1_strand(p->b)) {
                    plp_col->rv_counts[nt4] += 1;
               } else {
                    plp_col->fw_counts[nt4] += 1;
               }
                
               if (bi) {
                    char *t = (char*)(bi+1); /* 1 is type */
#if 0
                    int j;
                    printf("At %d qpos %d: BI=%s", plp_col->pos+1, p->qpos+1, t);
                    for (j = 0; j < p->indel; ++j) {
                         printf(" %c:%d-%d-%d", t[p->qpos+j], t[p->qpos+j-1]-33, t[p->qpos+j]-33, t[p->qpos+j+1]-33);
                    }
                    printf("\n");
#endif
                    /* adding 1 value representing whole insert */
                    PLP_COL_ADD_QUAL(& plp_col->ins_quals, t[p->qpos]);
               }

               if (bd) {
                    char *t = (char*)(bd+1);  /* 1 is type */
#if 0
                    int j;
                    printf("At %d qpos %d: BD=%sx", plp_col->pos+1, p->qpos+1, t);
                    for (j = 0; j < p->indel; ++j) {
                         printf(" %c:%d-%d-%d", t[p->qpos+j], t[p->qpos+j-1]-33, t[p->qpos+j]-33, t[p->qpos+j+1]-33);
                    }
                    printf("\n");
#endif
                    PLP_COL_ADD_QUAL(& plp_col->del_quals, t[p->qpos]);
                    /* adding 1 value representing whole del */
               }
          } /* ! p->is_del */
          
#if 0
          /* FIXME so what happens if is_del?
           *
           * special case p->is_refskip is
           * possible see ('spliced alignment'?)
           *
           * Observation: if indel<0 == del then they are followed by
           * is_del's. shouldn't that rather happen for ins? 
           */
          if (p->is_del || p->indel) {
               fflush(stdout); fflush(stderr);
               LOG_FIXME("At %s:%d p->is_del=%d p->is_refskip=%d p->indel=%d\n", 
                         plp_col->target, plp_col->pos+1, p->is_del, p->is_refskip, p->indel);
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

          if (p->indel != 0) {
               /* insertion (+)
                */
               if (p->indel > 0) {
                    plp_col->num_ins += 1;
#ifdef FIXME
                    LOG_FIXME("%s\n", "FIXME:need-to-save-ins-allele-and-its-counts.above-is-just-a-total.use-list/array-instead-of-hash.dont-expect-many-alleles");
#endif
                    plp_col->sum_ins += p->indel;
#ifdef PRINT_INDEL
                    putchar('+'); printw(p->indel, stdout);
                    putchar(' ');
                    int j;
                    for (j = 1; j <= p->indel; ++j) {
                         int c = bam_nt16_rev_table[bam1_seqi(bam1_seq(p->b), p->qpos + j)];
                         putchar(bam1_strand(p->b)? tolower(c) : toupper(c));
                    }
                    printf("\n");
#endif 
               /* deletion (-)
                */
               } else if (p->indel < 0) {
                    plp_col->num_dels += 1;
#ifdef FIXME
                    LOG_FIXME("%s\n", "FIXME:need-to-save-del-allele-and-its-counts.above-is-just-a-total.use-list/array-instead-of-hash.dont-expect-many-alleles");
#endif
                    plp_col->sum_dels -= p->indel;
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

     plp_col->coverage -= num_skips;
     
     /* determine consensus from 'counts'. will never produce N on tie  */
     plp_col->cons_base = bam_nt4_rev_table[
          argmax_d(base_counts, NUM_NT4)];

     if (debug) {
          plp_col_debug_print(plp_col, stderr);
     }
#if 0
     plp_col_mpileup_print(plp_col, conf, stdout);
#endif

     for (i = 0; i < NUM_NT4; ++i) {
          assert(plp_col->fw_counts[i] + plp_col->rv_counts[i] == plp_col->base_quals[i].n);
          assert(plp_col->base_quals[i].n == plp_col->map_quals[i].n);
#ifdef USE_SOURCEQUAL
           assert(plp_col->map_quals[i].n == plp_col->source_quals[i].n);
#endif
     }
}
/* compile_plp_col() */



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
    bam_header_t *h = 0;
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
     * files, whereas we allow only one. * however, if we keep it
     * close to the original source then a diff * against future
     * versions of samtools is easie
     *
     */
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
        data[i]->conf = mplp_conf;
        h_tmp = bam_header_read(data[i]->fp);
        if ( !h_tmp ) {
             fprintf(stderr,"[%s] fail to read the header of %s\n", __func__, fn[i]);
             exit(1);
        }
        data[i]->h = i? h : h_tmp; /* for i==0, "h" has not been set yet */

        if (mplp_conf->reg) {
            int beg, end;
            bam_index_t *idx;
            idx = bam_index_load(fn[i]);
            if (idx == 0) {
                fprintf(stderr, "[%s] fail to load index for %d-th input.\n", __func__, i+1);
                exit(1);
            }
            if (bam_parse_region(h_tmp, mplp_conf->reg, &tid, &beg, &end) < 0) {
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
    LOG_DEBUG("%s\n", "BAM header initialized");

    
    if (tid0 >= 0 && mplp_conf->fai) { /* region is set */
        ref = faidx_fetch_seq(mplp_conf->fai, h->target_name[tid0], 0, 0x7fffffff, &ref_len);
        LOG_DEBUG("%s\n", "sequence fetched");
        ref_tid = tid0;
        for (i = 0; i < n; ++i) data[i]->ref = ref, data[i]->ref_id = tid0;
    } else {
         ref_tid = -1;
         ref = 0;
    }
    iter = bam_mplp_init(n, mplp_func, (void**)data);
    max_depth = mplp_conf->max_depth;
    if (max_depth * 1 > 1<<20)
        fprintf(stderr, "(%s) Max depth is above 1M. Potential memory hog!\n", __func__);
    if (max_depth * 1 < 8000) {
        max_depth = 8000 / 1;
        fprintf(stderr, "<%s> Set max per-file depth to %d\n", __func__, max_depth);
    }
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
                 if (NULL == ref) {
                      LOG_FATAL("%s\n", "Given reference fasta file doesn't seem to contain the right sequence(s).");
                      return -1;
                 }
            }
            for (i = 0; i < n; ++i) 
                 data[i]->ref = ref, data[i]->ref_id = tid;
            ref_tid = tid;
        }
        i=0; /* i is 1 for first pos which is a bug due to the removal
              * of one of the loops, so reset here */

        plp_counter += 1;
        if (0 == plp_counter%100000) {
             LOG_VERBOSE("Still alive and happily crunching away on pos"
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
    bam_header_destroy(h);
    for (i = 0; i < n; ++i) {
        bam_close(data[i]->fp);
        if (data[i]->iter) bam_iter_destroy(data[i]->iter);
        free(data[i]);
    }
    free(data); free(plp); free(ref); free(n_plp);
    return 0;
}
/* mpileup() */
