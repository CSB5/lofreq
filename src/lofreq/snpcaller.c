/* -*- c-file-style: "k&r"; indent-tabs-mode: nil; -*-
 *
 * FIXME Copyright update
 *
 *********************************************************************/


#define TIMING 0

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include <ctype.h>
#include <float.h>

#include "fet.h"
#include "utils.h"
#include "log.h"

#include "snpcaller.h"

#if TIMING
#include <time.h>
#endif


/* median number of best hits in bwa for one WGS sample: 3. for one
 * exome sample it's 2. conservative choice is 3. FIXME make user
 * choice 
*/
#define MQ0_ERRPROB 0.66

#define LOGZERO -1e100 
/* FIXME shouldn't we use something from float.h ? */


#if 0
#define DEBUG
#endif

#if 0
#define TRACE
#endif

#if 0
#define NAIVE
#endif

/* scale 0-60 to from 0-254 and shrink
 * Y = 254.0/60.0 * MQ * (MQ**X)/(60**X)
 *
 * if 20 should stay 20
 * 20 = 254/60.0 * 20 * (20**X)/(60**X) 
 * 60/254.0 = (20**X)/(60**X)
 * (20/60.0)**X = 60/254.0
 * since a**x = y equals log_a(y) = x
 * x = log_a(60/254.0); a=20/60.0;
 * x = 1.3134658329243962 
 */
#if 0
#define SCALE_MQ 1
#endif
#define SCALE_MQ_FAC  1.3134658329243962

#if 0
#define TRUE_MQ_BWA_HG19_EXOME_2X100_SIMUL 1
#endif

/* filled in missing values with the min of the two neighbouring values */
static int MQ_TRANS_TABLE[61] = {
0, /* actually 1 but 0 should stay 0 */
1,
3,
4,
5,
5,
8,
9,
4,
8,
14,
17,
22,
25,
25,
29,
32,
33,
34,
34, /* NA */
34,
34, /* NA */
34, /* NA */
34,
34, /* NA */
34, /* NA */
34, /* NA */
34, /* NA */
34, /* NA */
41,
41, /* NA */
41, /* NA */
41, /* NA */
41, /* NA */
41, /* NA */
41, /* NA */
41, /* NA */
41, /* NA */
50,
46, /* NA */
46, /* NA */
46, /* NA */
46, /* NA */
46, /* NA */
46, /* NA */
46, /* NA */
46, /* NA */
46, /* NA */
46, /* NA */
46, /* NA */
46, /* NA */
46,
46, /* NA */
46, /* NA */
54,
37,
37, /* NA */
45,
45, /* NA */
45, /* NA */
67};



double log_sum(double log_a, double log_b);
double log_diff(double log_a, double log_b);
double probvec_tailsum(const double *probvec, int tail_startindex,
                       int probvec_len);
double *naive_calc_prob_dist(const double *err_probs, int N, int K);
double *pruned_calc_prob_dist(const double *err_probs, int N, int K, 
                      long long int bonf_factor, double sig_level);


static int mq_trans_range_violation_printed = 0;

int mq_trans(int mq) {
     if (mq<=60 && mq>=0) {
          return MQ_TRANS_TABLE[mq]; 

     } else if (mq!=255) {
          if (! mq_trans_range_violation_printed) {
               LOG_WARN("MQ value %d is outside of valid range defined in translation table\n", mq);
               mq_trans_range_violation_printed = 1;
          }
     }
     return mq;
}

/* J = PM  +  (1-PM) * PS  +  (1-PM) * (1-PS) * PA + (1-PM) * (1-PS) * (1-PA) * PB, where
   PJ = joined error probability
   PM = mapping prob.
   PS = source/genome err.prob.
   PS = mapping error profile prop.
   PB = base err.prob.
   
   Or in simple English:
   either this is a mapping error
   or
   not a mapping error, but a genome/source error
   or
   not mapping error and not genome/source error but a aligner error
   or
   not mapping error and not genome/source error and not aligner error but a base error
   
 * NOTE: for mq the standard says that 255 means NA. In this function we
 * use -1 instead, and treat 255 as valid phred-score so you might want
 * to change mq before calling this functio
 *
 */
double
merge_srcq_baseq_mapq_and_alnq(const int sq, const int bq, const int mq, const int aq)
{
     double sp, bp, mp, ap, jp; /* corresponding probs */
     

     bp = PHREDQUAL_TO_PROB(bq);

     if (-1 == mq) {
          mp = 0.0;
     } else if (0 == mq) {
          mp = MQ0_ERRPROB;
     } else {
          mp = PHREDQUAL_TO_PROB(mq);
     }

     if (-1 == sq) {
          sp = 0.0;
     } else {
          sp = PHREDQUAL_TO_PROB(sq);
     }

     if (-1 == aq) {
          ap = 0.0;
     } else {
          ap = PHREDQUAL_TO_PROB(aq);
     }

     jp = mp + (1.0 - mp) * sp + (1.0-mp) * (1.0-sp) * ap +  (1.0-mp) * (1.0-sp) * (1.0-ap) * bp;
#ifdef DEBUG
     LOG_DEBUG("jp=%g  =  mp=%g  +  (1.0-mp=%g)*sp=%g  +  (1-mp=%g)*(1-sp=%g)*ap=%g  + (1-mp=%g)*(1-sp=%g)*(1-ap=%g)*bp=%g\n",
               jp, mp, mp, sp, mp, sp, ap, mp, sp, ap, bp);
#endif

     return jp;
}



/* J = PM  +  (1-PM) * PS  +  (1-PM) * (1-PS) * PB, where
   PJ = joined error probability
   PM = mapping err.prob.
   PS = source/genome err.prob.
   PB = base err.prob.
   
   Or in simple English:
   either this is a mapping error
   or
   not a mapping error, but a genome/source error
   or
   not mapping error and no genome/source error AND base-error

 * NOTE: for mq the standard says that 255 means NA. In this function we
 * use -1 instead, and treat 255 as valid phred-score so you might want
 * to change mq before calling this functio
 *
*/
double
merge_srcq_baseq_and_mapq(const int sq, const int bq, const int mq)
{
     double sp, mp, bp, jp; /* corresponding probs */
     
     if (-1 == sq) {
          sp = 0.0;
     } else {
          sp = PHREDQUAL_TO_PROB(sq);
     }

     bp = PHREDQUAL_TO_PROB(bq);

     if (-1 == mq) {
          mp = 0.0;
     } else if (0 == mq) {
          mp = MQ0_ERRPROB;
     } else {
          mp = PHREDQUAL_TO_PROB(mq);
     }

     jp = mp + (1.0 - mp) * sp + (1-mp) * (1-sp) * bp;
#ifdef DEBUG
     LOG_DEBUG("jp = %g  =  mp=%g  +  (1.0 - mp=%g) * sp=%g  +  (1-mp=%g) * (1-sp=%g) * bp=%g\n", jp, mp, mp, sp, mp, sp, bp);
#endif

     return jp;
}


/* "Merge" MQ and BQ if requested using the following equation:
 *  P_jq = P_mq * + (1-P_mq) P_bq.
 *
 * NOTE: for mq the standard says that 255 means NA. In this function we
 * use -1 instead, and treat 255 as valid phred-score so you might want
 * to change mq before calling this functio
 *

 */
double
merge_baseq_and_mapq(const int bq, const int mq)
{
     double mp, bp, jp; /* corresponding probs */

     bp = PHREDQUAL_TO_PROB(bq);
     if (-1 == mq) {
          return bp;
     }

     if (0 == mq) {
          mp = MQ0_ERRPROB;
     } else {
          mp = PHREDQUAL_TO_PROB(mq);
     }

     /* No need to do computation in phred-space as
      * numbers won't get small enough.
      */

     /* note: merging Q1 with anything else will result in Q0. */
     jp = mp + (1.0 - mp) * bp;
#ifdef DEBUG
     LOG_DEBUG("P_M + (1-P_M) P_B:   %g + (1.0 - %g) * %g = %g  ==  Q%d + (1.0 - Q%d) * Q%d  =  Q%d\n",
               mp, mp, bp, jp, mq, mq, bq, PROB_TO_PHREDQUAL_SAFE(jp));
#endif
#if 0
     LOG_DEBUG("BQ %d after merging with MQ %d = %d\n", bq, mq, PROB_TO_PHREDQUAL_SAFE(jp));
#endif

     return jp;
}
/* merge_baseq_and_mapq() */



void
plp_to_errprobs(double **err_probs, int *num_err_probs, 
                int *alt_bases, int *alt_counts, int *alt_raw_counts,
                const plp_col_t *p, snvcall_conf_t *conf)
{
     int i, j;
     int alt_idx;
     int avg_qual = -1;

     if (NULL == ((*err_probs) = malloc(p->coverage * sizeof(double)))) {
          /* coverage = base-count after read level filtering */
          fprintf(stderr, "FATAL: couldn't allocate memory at %s:%s():%d\n",
                  __FILE__, __FUNCTION__, __LINE__);
          return;
     }

     /* determine avg_qual if needed, which is to be the median
      * quality of reference bases (ie. cons bases since vars are
      * always called against cons_base). compute average of all
      * unfiltered cons bases in this column (pretending we did that
      * many experiments)
     */
     if (-1 == conf->def_altbq) {

          for (i=0; i<NUM_NT4; i++) {
               int nt = bam_nt4_rev_table[i];
               if (nt != p->cons_base) {
                    continue;
               }

               if (p->base_quals[i].n == 0) {
                    /* set avg_qual to non-sense value, won't use it anyway
                     * if there were no 'errors'. can't return here though
                     * since we still want to report the consvar */
                    avg_qual = -1;
               } else {
                    int *ref_quals = malloc(sizeof(int) * p->base_quals[i].n);
                    memcpy(ref_quals, p->base_quals[i].data, sizeof(int) * p->base_quals[i].n);
                    avg_qual = int_median(ref_quals, p->base_quals[i].n);
                    free(ref_quals);
                    break; /* there can only be one */
               }
          }

     }


     (*num_err_probs) = 0;
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
               int bq = p->base_quals[i].data[j];                
               int mq = -1;
               int sq = -1;
#ifdef USE_ALNERRPROF
               int aq = -1;
#endif
               double final_err_prob; /* == final quality used for snv calling */
               
               if ((conf->flag & SNVCALL_USE_MQ) && p->map_quals[i].n) {
                    mq = p->map_quals[i].data[j];
                    /*according to spec 255 is unknown */
                    if (mq == 255) {
                         mq = -1;
                    }
#ifdef SCALE_MQ
                    mq = 254/60.0*mq * pow(mq, SCALE_MQ_FAC)/pow(60, SCALE_MQ_FAC);
#elif defined(TRUE_MQ_BWA_HG19_EXOME_2X100_SIMUL)
                    mq = mq_trans(mq); 
#endif
               }

#ifdef USE_SOURCEQUAL
               if (p->source_quals[i].n) {
                    sq = p->source_quals[i].data[j];
                    if (! (conf->flag & SNVCALL_USE_SQ)) {
                         sq = -1; /* i.e. NA */
                    }
               }
#endif

#ifdef USE_ALNERRPROF
               if (p->alnerr_qual[i].n) {
                    aq = p->alnerr_qual[i].data[j];
               }
#endif

#ifdef USE_ALNERRPROF
               final_err_prob = merge_srcq_baseq_mapq_and_alnq(sq, bq, mq, aq);
#else
               final_err_prob = merge_srcq_baseq_and_mapq(sq, bq, mq);
#endif

               /* special treatment of alt bases */
               if (is_alt_base) {
                    alt_raw_counts[alt_idx] += 1;

#if 0
                    LOG_FIXME("alt_base %d: bq=%d merged q=%d p=%f\n", alt_idx, bq, PROB_TO_PHREDQUAL_SAFE(final_err_prob), final_err_prob);
#endif
                    /* ignore if below bq or merged threshold */
                    if (bq < conf->min_altbq || PROB_TO_PHREDQUAL_SAFE(final_err_prob) < DEFAULT_MIN_ALT_MERGEDQ) {
                         continue; 
                    }
                    alt_counts[alt_idx] += 1;
                    /* if passed filter, set to default */
                    if (-1 == conf->def_altbq) {
                         /* ...change bq which also requires change of final_err_prob */
                         bq = avg_qual;
#ifdef USE_ALNERRPROF
                         final_err_prob = merge_srcq_baseq_mapq_and_alnq(sq, bq, mq, aq);
#else
                         final_err_prob = merge_srcq_baseq_and_mapq(sq, bq, mq);
#endif

                    } else {
                         /* easy case: just set to default */
                         final_err_prob = PHREDQUAL_TO_PROB(conf->def_altbq);
                    }
               }

#if 0
               LOG_FIXME("%s:%d %c bq=%d mq=%d finalq=%d is_alt_base=%d\n", p->target, p->pos+1, nt, bq, mq, PROB_TO_PHREDQUAL_SAFE(final_err_prob), is_alt_base);
#endif

               (*err_probs)[(*num_err_probs)++] = final_err_prob;
          }
     }
}



/* initialize members of preallocated snvcall_conf */
void
init_snvcall_conf(snvcall_conf_t *c) 
{
     memset(c, 0, sizeof(snvcall_conf_t));
     c->min_altbq = DEFAULT_MIN_ALT_BQ;
     c->def_altbq = DEFAULT_DEF_ALT_BQ; /* c->min_altbq; */
     c->min_cov = DEFAULT_MIN_COV;
     c->dont_skip_n = 0;
     c->bonf_dynamic = 1;
     c->bonf = 1;
     c->sig = 0.05;
     /* c->out = ; */
     c->flag = SNVCALL_USE_MQ;
}


void
dump_snvcall_conf(const snvcall_conf_t *c, FILE *stream) 
{
     fprintf(stream, "snvcall options\n");
     fprintf(stream, "  min_altbq      = %d\n", c->min_altbq);
     fprintf(stream, "  def_altbq      = %d\n", c->def_altbq);
     fprintf(stream, "  min_cov        = %d\n", c->min_cov);
     fprintf(stream, "  dont_skip_n    = %d\n", c->dont_skip_n);
     fprintf(stream, "  bonf           = %lld  (might get recalculated)\n", c->bonf);
     fprintf(stream, "  bonf_dynamic   = %d\n", c->bonf_dynamic);
     fprintf(stream, "  sig            = %f\n", c->sig);
/*     fprintf(stream, "  out            = %p\n", (void*)c->out);*/
     fprintf(stream, "  flag & SNVCALL_USE_MQ      = %d\n", c->flag&SNVCALL_USE_MQ?1:0);
     fprintf(stream, "  flag & SNVCALL_USE_SQ      = %d\n", c->flag&SNVCALL_USE_SQ?1:0);
     fprintf(stream, "  flag & SNVCALL_CONS_AS_REF = %d\n", c->flag&SNVCALL_CONS_AS_REF?1:0);
#ifdef SCALE_MQ
     LOG_WARN("%s\n", "MQ scaling switched on!");
#elif defined TRUE_MQ_BWA_HG19_EXOME_2X100_SIMUL
     LOG_WARN("%s\n", "MQ translation switched on!");
#endif
}



/**
 * @brief Computes log(exp(log_a) + exp(log_b))
 *
 * Taken from util.h of FAST source code:
 * http://www.cs.cornell.edu/~keich/FAST/fast.tar.gz
 * and using log1p
 */
double
log_sum(double log_a, double log_b)
{
    if (log_a > log_b) {
        return log_a + log1p(exp(log_b-log_a));
    } else {
        return log_b + log1p(exp(log_a-log_b));
    }
}
/* log_sum() */


/**
 * @brief Computes log(exp(log_a) - exp(log_b))
 *
 * Adapted from log_sum above and scala/breeze/numerics logDiff
 * See also http://stackoverflow.com/questions/778047/we-know-log-add-but-how-to-do-log-subtract
 *
 */
double
log_diff(double log_a, double log_b)
{
    if (log_a >= log_b) {
        return log_a + log1p(- exp(log_b-log_a));
    } else {
        return log_b + log1p(- exp(log_a-log_b));
    }
}
/* log_diff() */



/**
 * @brief Computes sum of probvec values (log space) starting from (including)
 * tail_startindex to (excluding) probvec_len
 *
 */
double
probvec_tailsum(const double *probvec, int tail_startindex, int probvec_len)
{
    double tailsum;
    int i;

    tailsum = probvec[tail_startindex];
    for (i=tail_startindex+1; i<probvec_len; i++) {
        tailsum = log_sum(tailsum, probvec[i]);
    }

    return tailsum;
}
/* probvec_tailsum() */


/**
 *
 */
double *
naive_calc_prob_dist(const double *err_probs, int N, int K)
{
     double *probvec = NULL;
     double *probvec_prev = NULL;
     double *probvec_swp = NULL;
     
     int n;
     fprintf(stderr, "CRITICAL(%s:%s:%d): Possibly buggy code. Use pruned_calc_prob_dist instead of me\n", 
             __FILE__, __FUNCTION__, __LINE__);
     exit(1);

    if (NULL == (probvec = malloc((N+1) * sizeof(double)))) {
        fprintf(stderr, "FATAL: couldn't allocate memory at %s:%s():%d\n",
                __FILE__, __FUNCTION__, __LINE__);
        return NULL;
    }
    if (NULL == (probvec_prev = malloc((N+1) * sizeof(double)))) {
        fprintf(stderr, "FATAL: couldn't allocate memory at %s:%s():%d\n",
                __FILE__, __FUNCTION__, __LINE__);
        free(probvec);
        return NULL;
    }

    /* init */
    probvec_prev[0] = 0.0; /* 0.0 = log(1.0) */

    for (n=1; n<N+1; n++) {
        int k;
        double log_pn, log_1_pn;
        double pn = err_probs[n-1];

        
        /* if pn=0 log(on) will fail. likewise if pn=1 (Q0) then
         * log1p(-pn) = log(1-1) = log(0) will fail. therefore test */
        if (fabs(pn) < DBL_EPSILON) {             
             log_pn = log(DBL_EPSILON);
        } else {
             log_pn = log(pn);
        }
        if (fabs(pn-1.0) < DBL_EPSILON) {             
             log_1_pn = log1p(-pn+DBL_EPSILON);
        } else {
             log_1_pn = log1p(-pn);
        }

#if 0
        fprintf(stderr, "DEBUG(%s:%s:%d): pn=%g log_pn=%g log_1_pn=%g err_probs[n=%d-1]=%g\n",
                __FILE__, __FUNCTION__, __LINE__, pn, log_pn, log_1_pn, n, err_probs[n-1]);
#endif

        k = 0;
        probvec[k] = probvec_prev[k] + log_1_pn;

        for (k=1; k<K; k++) {
             /* FIXME clang says: The left operand of '+' is a garbage value */
            probvec[k] = log_sum(probvec_prev[k] + log_1_pn,
                                 probvec_prev[k-1] + log_pn);
        }
        k = n;
        probvec[k] = probvec_prev[k-1] + log_pn;


        /* swap */
        probvec_swp = probvec;
        probvec = probvec_prev;
        probvec_prev = probvec_swp;
    }


    free(probvec_prev);    
    return probvec;
}
/* naive_prob_dist */



/**
 * Should really get rid of bonf_factor and sig_level here and
 * upstream as well
 *
 */
double *
pruned_calc_prob_dist(const double *err_probs, int N, int K, 
                      long long int bonf_factor, double sig_level)
{
    double *probvec = NULL;
    double *probvec_prev = NULL;
    double *probvec_swp = NULL;
    int n;

    if (NULL == (probvec = malloc((K+1) * sizeof(double)))) {
        fprintf(stderr, "FATAL: couldn't allocate memory at %s:%s():%d\n",
                __FILE__, __FUNCTION__, __LINE__);
        return NULL;
    }
    if (NULL == (probvec_prev = malloc((K+1) * sizeof(double)))) {
        fprintf(stderr, "FATAL: couldn't allocate memory at %s:%s():%d\n",
                __FILE__, __FUNCTION__, __LINE__);
        free(probvec);
        return NULL;
    }

    for (n=0; n<N; n++) {
         /*LOG_FIXME("err_probs[n=%d]=%g\n", n, err_probs[n]);*/
         assert(err_probs[n] + DBL_EPSILON >= 0.0 && err_probs[n] - DBL_EPSILON <= 1.0);
    }

#ifdef DEBUG
    for (n=0; n<K+1; n++) {
        probvec_prev[n] = probvec[n] = 666.666;
    }
#endif

    /* init */
    probvec_prev[0] = 0.0; /* log(1.0) */

    for (n=1; n<=N; n++) {
        int k;
        double pvalue;
        double pn = err_probs[n-1];
        double log_pn, log_1_pn;


        /* if pn=0 log(on) will fail. likewise if pn=1 (Q0) then
         * log1p(-pn) = log(1-1) = log(0) will fail. therefore test */
        if (fabs(pn) < DBL_EPSILON) {             
             log_pn = log(DBL_EPSILON);
        } else {
             log_pn = log(pn);
        }
        if (fabs(pn-1.0) < DBL_EPSILON) {             
             log_1_pn = log1p(-pn+DBL_EPSILON);
        } else {
             log_1_pn = log1p(-pn);/* 0.0 = log(1.0) */
        }
        
#ifdef TRACE
		fprintf(stderr, "DEBUG(%s:%s:%d): n=%d err_probs[n-1]=%g pn=%g log_pn=%g log_1_pn=%g\n", 
                __FILE__, __FUNCTION__, __LINE__, n, err_probs[n-1], pn, log_pn, log_1_pn);
#endif

        if(n < K) {
            probvec_prev[n] = LOGZERO;
        }

        for (k=MIN(n,K-1); k>=1; k--) {
            assert(probvec_prev[k]<=0.0 && probvec_prev[k-1]<=0.0);
            probvec[k] = log_sum(probvec_prev[k] + log_1_pn,
                                 probvec_prev[k-1] + log_pn);            
        }
        k = 0;
        assert(probvec_prev[k]<=0.0);
        probvec[k] = probvec_prev[k] + log_1_pn;

#ifdef TRACE
        for (k=0; k<=MIN(n, K-1); k++) {
            fprintf(stderr, "DEBUG(%s:%s:%d): probvec[k=%d] = %g\n", 
                    __FILE__, __FUNCTION__, __LINE__, k, probvec[k]);
        }
        for (k=0; k<=MIN(n,K-1); k++) {
            fprintf(stderr, "DEBUG(%s:%s:%d): probvec_prev[k=%d] = %g\n", 
                    __FILE__, __FUNCTION__, __LINE__, k, probvec_prev[k]);
        }
#endif

        if (n==K) {
            probvec[K] = probvec_prev[K-1] + log_pn;
            /* FIXME check here as well */

        } else if (n > K) { 
             /*LOG_FIXME("probvec_prev[K=%d]=%g probvec_prev[K=%d -1]=%g\n", K, probvec_prev[K], K, probvec_prev[K-1]);*/
             assert(probvec_prev[K]-DBL_EPSILON<=0.0 && probvec_prev[K-1]-DBL_EPSILON<=0.0);
             probvec[K] = log_sum(probvec_prev[K], probvec_prev[K-1]+log_pn);
             pvalue = exp(probvec[K]);
             
             if (pvalue * (double)bonf_factor >= sig_level) {
#ifdef DEBUG
                  fprintf(stderr, "DEBUG(%s:%s:%d): early exit at n=%d with pvalue %g\n", 
                          __FILE__, __FUNCTION__, __LINE__, n, pvalue);
#endif
                  goto free_and_exit;
             }
        }

        assert(! isinf(probvec[0])); /* used to happen when first q=0 */

        /* swap */
        probvec_swp = probvec;
        probvec = probvec_prev;
        probvec_prev = probvec_swp;
    }

 free_and_exit:
    free(probvec_prev);    

    return probvec;
}
/* pruned_calc_prob_dist */


#ifdef PSEUDO_BINOMIAL
/* binomial test using poissbin. only good for high n and small prob.
 * returns -1 on error */
int
pseudo_binomial(double *pvalue, 
                int num_success, int num_trials, double succ_prob)
{
     const long long int bonf = 1.0;
     const double sig = 1.0;
     double *probvec = NULL;
     double *probs;
     int i;

     fprintf(stderr, "WARNING(%s): this function only approximates the binomial for high n and small p\n", __FUNCTION__);
     if (NULL == (probs = malloc((num_trials) * sizeof(double)))) {
          fprintf(stderr, "FATAL: couldn't allocate memory at %s:%s():%d\n",
                  __FILE__, __FUNCTION__, __LINE__);
          return -1;
     }

     for (i=0; i<num_trials; i++) {
          probs[i] = succ_prob;
     }     

     probvec = poissbin(pvalue, probs,
                        num_trials, num_success,
                        bonf, sig);
     free(probvec);
     free(probs);

     return 0;
}
#endif



/* main logic. return of probvec (needs to be freed by caller allows
 * to check pvalues for other numbers < (original num_failures), like
 * so: exp(probvec_tailsum(probvec, smaller_numl, orig_num+1)) but
 * only if first pvalue was below limits implied by bonf and sig.
 * default pvalue is DBL_MAX (1 might still be significant).
 *
 *  note: pvalues > sig/bonf are not computed properly
 */       
double *
poissbin(double *pvalue, const double *err_probs,
         const int num_err_probs, const int num_failures, 
         const long long int bonf, const double sig) 
{
    double *probvec = NULL;
#if TIMING
    clock_t start = clock();
    int msec;
#endif
    *pvalue = DBL_MAX;

#if TIMING
    start = clock();
#endif
#ifdef NAIVE
    probvec = naive_prob_dist(err_probs, num_err_probs,
                                    num_failures);
#else
    probvec = pruned_calc_prob_dist(err_probs, num_err_probs,
                                    num_failures, bonf, sig);    
#endif
#if TIMING
    msec = (clock() - start) * 1000 / CLOCKS_PER_SEC;
    fprintf(stderr, "calc_prob_dist() took %d s %d ms\n", msec/1000, msec%1000);
#endif

    *pvalue = exp(probvec[num_failures]); /* no need for tailsum here */
    assert(! isnan(*pvalue));
    return probvec;
}



/**
 * @brief
 * 
 * pvalues computed for each of the NUM_NONCONS_BASES noncons_counts
 * will be written to snp_pvalues in the same order. If pvalue was not
 * computed (always insignificant) its value will be set to DBL_MAX
 * 
 */
int
snpcaller(double *snp_pvalues, 
          const double *err_probs, const int num_err_probs, 
          const int *noncons_counts, 
          const long long int bonf_factor, const double sig_level)
{
    double *probvec = NULL;
    int i;
    int max_noncons_count = 0;
    double pvalue;

#ifdef DEBUG
    fprintf(stderr, "DEBUG(%s:%s():%d): num_err_probs=%d noncons_counts=%d,%d,%d bonf_factor=%lld sig_level=%f\n", 
            __FILE__, __FUNCTION__, __LINE__, 
            num_err_probs, noncons_counts[0], noncons_counts[1], noncons_counts[2],
            bonf_factor, sig_level);
#endif

    /* initialise empty results so that we can return anytime */
    for (i=0; i<NUM_NONCONS_BASES; i++) {
        snp_pvalues[i] = DBL_MAX;
    }
    
    /* determine max non-consensus count */
    for (i=0; i<NUM_NONCONS_BASES; i++) {
        if (noncons_counts[i] > max_noncons_count) {
            max_noncons_count = noncons_counts[i];
        }
    }

    /* no need to do anything if no snp bases */
    if (0==max_noncons_count) {
        goto free_and_exit;
    }

    probvec = poissbin(&pvalue, err_probs, num_err_probs,
                       max_noncons_count, bonf_factor, sig_level);
    
    if (pvalue * (double)bonf_factor >= sig_level) {
#ifdef DEBUG
        fprintf(stderr, "DEBUG(%s:%s():%d): Most frequent SNV candidate already gets not signifcant pvalue of %g * %lld >= %g\n", 
                __FILE__, __FUNCTION__, __LINE__, 
                pvalue, bonf_factor, sig_level);
#endif
        goto free_and_exit;
    }


    /* report p-value for each non-consensus base
     */
#if 0
    for (i=1; i<max_noncons_count+1; i++) {        
        fprintf(stderr, "DEBUG(%s:%s():%d): prob for count %d=%g\n", 
                __FILE__, __FUNCTION__, __LINE__, 
                i, exp(probvec[i]));
    }
#endif
#if 0
    for (i=1; i<max_noncons_count+1; i++) {        
        fprintf(stderr, "DEBUG(%s:%s():%d): pvalue=%g for noncons_counts %d\n", 
                __FILE__, __FUNCTION__, __LINE__, 
                exp(probvec_tailsum(probvec, i, max_noncons_count+1)), i);
    }
#endif

    for (i=0; i<NUM_NONCONS_BASES; i++) { 
        if (0 != noncons_counts[i]) {
            pvalue = exp(probvec_tailsum(probvec, noncons_counts[i], max_noncons_count+1));
            snp_pvalues[i] = pvalue;
#ifdef DEBUG
            fprintf(stderr, "DEBUG(%s:%s():%d): i=%d noncons_counts=%d max_noncons_count=%d pvalue=%g\n", 
                    __FILE__, __FUNCTION__, __LINE__, 
                    i, noncons_counts[i], max_noncons_count, pvalue);                  
#endif
        }
    }

 free_and_exit:
    if (NULL != probvec) {
        free(probvec);
    }

    return 0;
}
/* snpcaller() */


#ifdef SNPCALLER_MAIN


/* 
 * gcc -pedantic -Wall -g -std=gnu99 -O2 -DSNPCALLER_MAIN -o snpcaller snpcaller.c utils.c log.c
 * 
 * to test the pvalues (remember: large n and small p when comparing to binomial)
 *
 * >>> [scipy.stats.binom_test(x, 10000, 0.0001) for x in [2, 3, 4]]
 * [0.26424111735042727, 0.080292199242652212, 0.018982025450177534]
 * 
 * ./snpcaller 4 10000 0.0001
 * prob from snpcaller(): (.. 0.264204 .. 0.0802738 ..) 0.0189759
 *
 * 
 * to test probs:

 * >>> print(zip(range(1,11), scipy.stats.binom.pmf(range(1,11), 10000, 0.0001)))
 * [(1, 0.36789783621841865), (2, 0.18394891810853542), (3, 0.061310173792593993), (4, 0.015324477633007457), (5, 0.0030639759659701437), (6, 0.00051045837549934915), (7, 7.2886160113568433e-05), (8, 9.1053030054241777e-06), (9, 1.0109920728573426e-06), (10, 1.0101831983287595e-07)]
 *
 * ./snpcaller 10 10000 0.0001
 * prob for count 1=...
 *
 */
int main(int argc, char *argv[]) {
     int num_success;
     int num_trials;
     double succ_prob;
     verbose = 1;

     if (argc<4) {
          LOG_ERROR("%s\n", "need num_success num_trials and succ_prob as args");
          return -1;
     }
     num_success = atoi(argv[1]);
     num_trials = atoi(argv[2]);
     succ_prob = atof(argv[3]);

     LOG_VERBOSE("num_success=%d num_trials=%d succ_prob=%f\n", num_success, num_trials, succ_prob);


#ifdef PSEUDO_BINOMIAL
     {
          double pvalue;
          if (-1 == pseudo_binomial(&pvalue, 
                                    num_success, num_trials, succ_prob)) {
               LOG_ERROR("%s\n", "pseudo_binomial() failed");
               return -1;
          }
          printf("pseudo_binomial: %g\n", pvalue);
     }
#endif


#if 1
     {
          double snp_pvalues[NUM_NONCONS_BASES];
          int noncons_counts[NUM_NONCONS_BASES];
          double *err_probs;
          int i;

          if (NULL == (err_probs = malloc((num_trials) * sizeof(double)))) {
               fprintf(stderr, "FATAL: couldn't allocate memory at %s:%s():%d\n",
                       __FILE__, __FUNCTION__, __LINE__);
               return -1;
          }
          for (i=0; i<num_trials; i++) {
               err_probs[i] = succ_prob;
          }

          noncons_counts[0] = num_success;
          noncons_counts[1] = num_success-1;
          noncons_counts[2] = num_success-2;

          snpcaller(snp_pvalues, err_probs, num_trials, noncons_counts, 1, 1);

          printf("prob from snpcaller(): (.. -2:%g .. -1:%g ..) %g\n", snp_pvalues[2], snp_pvalues[1], snp_pvalues[0]);
          free(err_probs);
     }
#endif
}
#endif
