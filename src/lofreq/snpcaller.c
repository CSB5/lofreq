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



#define TIMING 0

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include <ctype.h>
#include <float.h>
#include <errno.h>
#include <fenv.h>

#include "fet.h"
#include "utils.h"
#include "log.h"

#include "snpcaller.h"
#if TIMING
#include <time.h>
#endif


/* Converting MQ=0 into prob would 'kill' a read. Previously used 0.66 here since
   the median number of best hits in BWA for one examined human wgs sample
   was 3 (sadly BWA-MEM doesn't produce X0 tags anymore). For simplicity's
   sake, give MQ0 the benefit of doubt and assume that one only one other best
   location existed, i.e. use 0.5
*/
#define MQ0_ERRPROB 0.5

#define LOGZERO -1e100
/* shouldn't we use something from float.h ? */


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
#define SCALE_MQ_FAC  1.3134658329243962
#endif

#if 0
/* filled in missing values with the min of the two neighbouring values */
static int MQ_TRANS_BWA_062_SAMPE_HG19_2X100_SIMUL[61] = {
1,
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


static int MQ_TRANS_BWA_079_MEM_HG19_CHR22_2X75_SIMUL[72] = {
     19,
     48,
     66,
     47,
     58,
     59,
     58,
     55,
     60,
     54,
     57,
     58,
     58,
     58,
     65,
     56,
     60,
     62,
     62,
     57,
     57,
     54,
     50,
     51,
     52,
     49,
     74,
     49,
     49,
     77,
     77,
     77,
     77,
     68,
     68,
     74,
     74,
     74,
     74,
     74,
     68,
     68,
     68,
     68,
     71,
     71,
     71,
     71,
     77,
     77,
     77,
     77,
     77,
     77,
     77,
     77,
     77,
     74,
     74,
     74,
     77,
     77,
     77,
     77,
     77,
     77,
     77,
     77,
     77,
     77,
     77};

#if 0
#define MQ_TRANS_TABLE MQ_TRANS_BWA_062_SAMPE_HG19_2X100_SIMUL
#else
#define MQ_TRANS_TABLE MQ_TRANS_BWA_079_MEM_HG19_CHR22_2X75_SIMUL
#endif
static int mq_trans_range_violation_printed = 0;

#endif


double log_sum(double log_a, double log_b);
double log_diff(double log_a, double log_b);
double probvec_tailsum(const double *probvec, int tail_startindex,
                       int probvec_len);
double *naive_calc_prob_dist(const double *err_probs, int N, int K);
double *pruned_calc_prob_dist(const double *err_probs, int N, int K,
                      long long int bonf_factor, double sig_level);



#ifdef MQ_TRANS_TABLE
int mq_trans(int mq) {
#if 1
     if (mq<=72  && mq>=0) {
#else
     if (mq<=60  && mq>=0) {
#endif
          return MQ_TRANS_TABLE[mq];

     } else if (mq!=255) {
          if (! mq_trans_range_violation_printed) {
               LOG_WARN("MQ value %d is outside of valid range defined in translation table\n", mq);
               mq_trans_range_violation_printed = 1;
          }
     }
     return mq;
}
#endif



/* PJ = PM + (1-PM)*PS + (1-PM)*(1-PS)*PA + (1-PM)*(1-PS)*(1-PA)*PB, where
 * PJ = joined error prob
 * PM = mapping error prob
 * PS = source/genome error prob
 * PA = base alignment error prob (BAQ)
 * PB = base error prob
 * Or in plain English:
 * either this is a mapping error
 * or
 * not, but a genome/source error
 * or
 * none of the above, but a base-alignment error
 * or
 * none of the above but a base-error
 *
 * In theory PS should go first but the rest is hard to compute then.
 * Using PM things get tractable and it intrinsically takes care of
 * PS.
 *
 * NOTE: the standard says that MQ=255 means NA. In this function we
 * use -1 instead for all unknown values, and treat 255 as valid
 * phred-score so you might want to change mq before.
 *
 */
double
merge_srcq_mapq_baq_and_bq(const int sq, const int mq, const int baq, const int bq)
{
     double sp, mp, bap, bp, jp; /* corresponding probs */

     if (-1 == sq) {
          sp = 0.0;
     } else {
          sp = PHREDQUAL_TO_PROB(sq);
     }

     if (-1 == mq) {
          mp = 0.0;
     } else if (0 == mq) {
          mp = MQ0_ERRPROB;
     } else {
          mp = PHREDQUAL_TO_PROB(mq);
     }

     if (-1 == baq) {
          bap = 0.0;
     } else {
          bap = PHREDQUAL_TO_PROB(baq);
     }

     if (-1 == bq) {
          bp = 0.0;
     } else {
          bp = PHREDQUAL_TO_PROB(bq);
     }

     /* FIXME do calculations in log space and return Q instead of p */
     jp = mp + (1.0-mp)*sp + (1-mp)*(1-sp)*bap + (1-mp)*(1-sp)*(1-bap)*bp;

#if 0
     LOG_DEBUG("sq=%d|%f mq=%d|%f baq=%d|%f bq=%d|%f. returning %f\n",
              sq, sp, mq, mp, baq, bap, bq, bp, jp);
#endif
     return jp;
}



void
plp_to_errprobs(double **err_probs, int *num_err_probs,
                int *alt_bases, int *alt_counts, int *alt_raw_counts,
                const plp_col_t *p, varcall_conf_t *conf)
{
     int i, j;
     int alt_idx;
     int avg_ref_bq = -1;

     if (NULL == ((*err_probs) = malloc(p->coverage_plp * sizeof(double)))) {
          /* coverage = base-count after read level filtering */
          fprintf(stderr, "FATAL: couldn't allocate memory at %s:%s():%d\n",
                  __FILE__, __FUNCTION__, __LINE__);
          return;
     }

     /* determine median ref bq in advance if needed
     */
     if (-1 == conf->def_alt_bq) {
          avg_ref_bq = -1;
          for (i=0; i<NUM_NT4; i++) {
               int nt = bam_nt4_rev_table[i];
               if (nt != p->ref_base) {
                    continue;
               }
               if (p->base_quals[i].n) {
                    int *ref_quals = malloc(sizeof(int) * p->base_quals[i].n);
                    memcpy(ref_quals, p->base_quals[i].data, sizeof(int) * p->base_quals[i].n);
                    avg_ref_bq = int_median(ref_quals, p->base_quals[i].n);
                    free(ref_quals);
                    break; /* there can only be one */
               }
          }
          LOG_DEBUG("avg_ref_bq=%d\n", avg_ref_bq);
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
          if (nt != p->ref_base) {
               is_alt_base = 1;
               alt_idx += 1;
               alt_bases[alt_idx] = nt;
               alt_counts[alt_idx] = 0;
               alt_raw_counts[alt_idx] = 0;
          }

          for (j=0; j<p->base_quals[i].n; j++) {
               int bq = -1;
               int mq = -1;
               int sq = -1;
               int baq = -1;
#ifdef USE_ALNERRPROF
               int aq = -1;
               LOG_FATAL("%s\n", "ALNERRPROF not supported anymore\n"); exit(1);
#endif
               double merged_err_prob; /* final quality used for snv calling */
               int merged_qual;

               if (p->base_quals[i].n) {
                    bq = p->base_quals[i].data[j];

                    /* bq filtering for all
                     * FIXME min_bq was meant to just influence quality calculations, not af etc
                     * but the filtering leaks out later, because alt_counts is computed and returned after filtering
                     */
                    if (bq < conf->min_bq) {
                         continue;
                    }

                    /* alt bq threshold and overwrite if needed */
                    if (is_alt_base) {
                         alt_raw_counts[alt_idx] += 1;
                         /* ignore altogether if below alt bq threshold */
                         if (bq < conf->min_alt_bq) {
                              continue;
                         } else if (-1 == conf->def_alt_bq)  {
                              bq = avg_ref_bq;
                         } else if (0 != conf->def_alt_bq)  {
                              bq = conf->def_alt_bq;
                         }
                         /* 0: keep original */
                    }
               }

               if ((conf->flag & VARCALL_USE_BAQ) && p->baq_quals[i].n) {
                    baq = p->baq_quals[i].data[j];
               }

               if ((conf->flag & VARCALL_USE_MQ) && p->map_quals[i].n) {
                    mq = p->map_quals[i].data[j];
                    /*according to spec 255 is unknown */
                    if (mq == 255) {
                         mq = -1;
                    }
#ifdef SCALE_MQ
                    mq = 254/60.0*mq * pow(mq, SCALE_MQ_FAC)/pow(60, SCALE_MQ_FAC);
#elif defined(MQ_TRANS_TABLE)
                    mq = mq_trans(mq);
#endif
               }

               if ((conf->flag & VARCALL_USE_SQ) && p->source_quals[i].n) {
                    sq = p->source_quals[i].data[j];
               }

               merged_err_prob = merge_srcq_mapq_baq_and_bq(sq, mq, baq, bq);
               merged_qual =  PROB_TO_PHREDQUAL_SAFE(merged_err_prob);

               /* min merged q filtering for all */
               if (merged_qual < conf->min_jq) {
                    continue;
               }

               if (is_alt_base) {
#if 0
                    LOG_debug("alt_base %d: bq=%d merged q=%d p=%f\n",
                              alt_idx, bq, PROB_TO_PHREDQUAL_SAFE(merged_err_prob), merged_err_prob);
#endif
                    /* apply alt merged qual threshold and overwrite if needed
                     */
                    if (merged_qual < conf->min_alt_jq) {
                         continue;
                    } else if (-1 == conf->def_alt_jq)  {
                         LOG_FATAL("%s\n", "median off ref joined q not implemented yet (FIXME)");
                         exit(1);
                    } else if (0 != conf->def_alt_jq)  {
                         merged_err_prob = PHREDQUAL_TO_PROB(conf->def_alt_jq);
                    }
                    /* 0: keep original */
                    alt_counts[alt_idx] += 1;
               }
               (*err_probs)[(*num_err_probs)++] = merged_err_prob;

#if 0
               LOG_FIXME("%s:%d %c bq=%d mq=%d finalq=%d is_alt_base=%d\n", p->target, p->pos+1, nt, bq, mq, PROB_TO_PHREDQUAL_SAFE(merged_err_prob), is_alt_base);
#endif
          }
     }
}

/* FIXME merge with plp_to_del_errprobs */
void
plp_to_ins_errprobs(double **err_probs, int *num_err_probs,
                    const plp_col_t *p, varcall_conf_t *conf,
                    char key[MAX_INDELSIZE]){

     if (NULL == ((*err_probs) = malloc(p->coverage_plp * sizeof(double)))) {
          /* coverage = base-count after read level filtering */
          fprintf(stderr, "FATAL: couldn't allocate memory at %s:%s():%d\n",
                  __FILE__, __FUNCTION__, __LINE__);
          free(err_probs);
          return;
     }

     (*num_err_probs) = 0;
     int i, j;
     double final_err_prob;
     int iq, aq, mq, sq;
     iq = aq = mq = sq = -1;

     for (i = 0; i < p->ins_quals.n; i++) {
          iq = mq = -1;
          iq = p->ins_quals.data[i];
          if (conf->flag & VARCALL_USE_MQ) {
               mq = p->ins_map_quals.data[i];
          }
          final_err_prob = merge_srcq_mapq_baq_and_bq(-1, mq, -1, iq);
          (*err_probs)[(*num_err_probs)++] = final_err_prob;
     }

     ins_event *it, *it_tmp;
     HASH_ITER(hh_ins, p->ins_event_counts, it, it_tmp) {
          for (j = 0; j < it->ins_quals.n; j++) {
               iq = aq = mq = sq = -1;
               iq = it->ins_quals.data[j];

               /* don't use idaq if not wanted or if not indel in question (FIXME does the latter amek sense)? */
               if ((conf->flag & VARCALL_USE_IDAQ) && (0 == strcmp(it->key, key))) {
                    aq = it->ins_aln_quals.data[j];
               }

               if ((conf->flag & VARCALL_USE_MQ) && it->ins_map_quals.n) {
                    mq = it->ins_map_quals.data[j];
                    /*according to spec 255 is unknown */
                    if (mq == 255) {
                         mq = -1;
                    }
               }

               if ((conf->flag & VARCALL_USE_SQ) && it->ins_source_quals.n)  {
                    sq = it->ins_source_quals.data[j];
               }
               
               final_err_prob = merge_srcq_mapq_baq_and_bq(sq, mq, aq, iq);
#ifdef TRACE
               LOG_DEBUG("+%s IQ:%d IAQ:%d MQ:%d SQ:%d EP:%lg\n",
                         it->key, iq, aq, mq, sq, final_err_prob);
#endif
               (*err_probs)[(*num_err_probs)++] = final_err_prob;
          }
     }
}

/* FIXME merge with plp_to_ins_errprobs */
void
plp_to_del_errprobs(double **err_probs, int *num_err_probs,
                    const plp_col_t *p, varcall_conf_t *conf,
                    char key[MAX_INDELSIZE]){
     if (NULL == ((*err_probs) = malloc(p->coverage_plp * sizeof(double)))) {
          /* coverage = base-count after read level filtering */
          fprintf(stderr, "FATAL: couldn't allocate memory at %s:%s():%d\n",
                  __FILE__, __FUNCTION__, __LINE__);
          free(err_probs);
          return;
     }

     (*num_err_probs) = 0;
     int i, j;
     double final_err_prob;
     int dq, aq, mq, sq;
     dq = aq = mq = sq = -1;

     for (i = 0; i < p->del_quals.n; i++) {
          dq = mq = -1;
          dq = p->del_quals.data[i];
          if (conf->flag & VARCALL_USE_MQ) {
               mq = p->del_map_quals.data[i];
          }
          final_err_prob = merge_srcq_mapq_baq_and_bq(-1, mq, -1, dq);
          (*err_probs)[(*num_err_probs)++] = final_err_prob;
     }

     del_event *it, *it_tmp;
     HASH_ITER(hh_del, p->del_event_counts, it, it_tmp) {
          for (j = 0; j < it->del_quals.n; j++) {
               dq = aq = mq = sq = -1;
               dq = it->del_quals.data[j];

               /* don't use idaq if not wanted or if not indel in question (FIXME does the latter amek sense)? */
               if ((conf->flag & VARCALL_USE_IDAQ) && (0 == strcmp(it->key, key))) {
                    aq = it->del_aln_quals.data[j];
               }

               if ((conf->flag & VARCALL_USE_MQ) && it->del_map_quals.n) {
                    mq = it->del_map_quals.data[j];
                    /*according to spec 255 is unknown */
                    if (mq == 255) {
                         mq = -1;
                    }
               }

               if ((conf->flag & VARCALL_USE_SQ) && it->del_source_quals.n) {
                         sq = it->del_source_quals.data[j];
               }

               final_err_prob = merge_srcq_mapq_baq_and_bq(sq, mq, aq, dq);
#ifdef TRACE
               LOG_DEBUG("+%s DQ:%d DAQ:%d MQ:%d SQ:%d EP:%lg\n",
                         it->key, dq, aq, mq, sq, final_err_prob);
#endif
               (*err_probs)[(*num_err_probs)++] = final_err_prob;
          }
     }
}

/* initialize members of preallocated varcall_conf */
void
init_varcall_conf(varcall_conf_t *c)
{
     memset(c, 0, sizeof(varcall_conf_t));

     c->min_bq = DEFAULT_MIN_BQ;
     c->min_alt_bq = DEFAULT_MIN_ALT_BQ;
     c->def_alt_bq = DEFAULT_DEF_ALT_BQ;

     c->min_jq = DEFAULT_MIN_JQ;
     c->min_alt_jq = DEFAULT_MIN_ALT_JQ;
     c->def_alt_jq = DEFAULT_DEF_ALT_JQ;

     c->min_cov = DEFAULT_MIN_COV;
     c->bonf_dynamic = 1;
     c->bonf_subst = 1;
     c->bonf_indel = 1;
     c->sig = DEFAULT_SIG;
     /* c->out = ; */
     c->flag |= VARCALL_USE_MQ;
     c->flag |= VARCALL_USE_BAQ;
     c->flag |= VARCALL_USE_IDAQ;
     c->only_indels = 0;
     c->no_indels = 0;
}


void
dump_varcall_conf(const varcall_conf_t *c, FILE *stream)
{
     fprintf(stream, "snvcall options\n");
     fprintf(stream, "  min_bq         = %d\n", c->min_bq);
     fprintf(stream, "  min_alt_bq     = %d\n", c->min_alt_bq);
     fprintf(stream, "  def_alt_bq     = %d\n", c->def_alt_bq);
     fprintf(stream, "  min_jq         = %d\n", c->min_jq);
     fprintf(stream, "  min_alt_jq     = %d\n", c->min_alt_jq);
     fprintf(stream, "  def_alt_jq     = %d\n", c->def_alt_jq);
     fprintf(stream, "  min_cov        = %d\n", c->min_cov);
     fprintf(stream, "  bonf_subst       = %lld  (might get recalculated)\n", c->bonf_subst);
     fprintf(stream, "  bonf_indel     = %lld  (might get recalculated)\n", c->bonf_indel);
     fprintf(stream, "  bonf_dynamic   = %d\n", c->bonf_dynamic);
     fprintf(stream, "  sig            = %f\n", c->sig);
/*     fprintf(stream, "  out            = %p\n", (void*)c->out);*/
     fprintf(stream, "  flag & VARCALL_USE_BAQ     = %d\n", c->flag&VARCALL_USE_BAQ?1:0);
     fprintf(stream, "  flag & VARCALL_USE_MQ      = %d\n", c->flag&VARCALL_USE_MQ?1:0);
     fprintf(stream, "  flag & VARCALL_USE_SQ      = %d\n", c->flag&VARCALL_USE_SQ?1:0);
     fprintf(stream, "  flag & VARCALL_USE_IDAQ    = %d\n", c->flag&VARCALL_USE_IDAQ?1:0);
#ifdef SCALE_MQ
     LOG_WARN("%s\n", "MQ scaling switched on!");
#elif defined MQ_TRANS_TABLE
     LOG_WARN("%s\n", "MQ translation switched on!");
#endif
     fprintf(stream, "  only_indels    = %d\n", c->only_indels);
     fprintf(stream, "  no_indels      = %d\n", c->no_indels);
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
             /* FIXME prune here as well? */

        } else if (n > K) {
             long double pvalue;
             int errsv = 0;
             /*LOG_FIXME("probvec_prev[K=%d]=%g probvec_prev[K=%d -1]=%g\n", K, probvec_prev[K], K, probvec_prev[K-1]);*/
             assert(probvec_prev[K]-DBL_EPSILON<=0.0 && probvec_prev[K-1]-DBL_EPSILON<=0.0);

             probvec[K] = log_sum(probvec_prev[K], probvec_prev[K-1]+log_pn);

             errno = 0;
             feclearexcept(FE_ALL_EXCEPT);

             pvalue = expl(probvec[K]);

             errsv = errno;
             if (errsv || fetestexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW | FE_UNDERFLOW)) {
                  if (pvalue < DBL_EPSILON) {
                       pvalue = LDBL_MIN;/* to zero but prevent actual 0 value */
                  } else {
                       pvalue = LDBL_MAX; /* might otherwise be set to 1 which might pass filters */
                  }
             }
             /* store as phred scores instead:

              Q = -10*log_10(e^X), where X=probvec[K]
              remember, log_b(x) = log_k(x)/log_k(b), i.e. log_10(Y) = log_e(Y)/log_e(10)
              therefore, Q = -10 * log_e(e^X)/log_e(10) = -10 * X/log_e(10)
              e.g.
              >>> from math import log, log10, e
              >>> X = -100
              >>> -10 * log10(e**X)
              434.29448190325184
              >>> -10 * X/log(10)
              434.2944819032518
             */
             if (pvalue * (double)bonf_factor > sig_level) {
#ifdef DEBUG
                  fprintf(stderr, "DEBUG(%s:%s:%d): early exit at n=%d K=%d with pvalue %Lg\n",
                          __FILE__, __FUNCTION__, __LINE__, n, K, pvalue);
#endif
                  free(probvec_prev);
                  return probvec;
             }
        }

        assert(! isinf(probvec[0])); /* used to happen when first q=0 */

        /* swap */
        probvec_swp = probvec;
        probvec = probvec_prev;
        probvec_prev = probvec_swp;
    }

    /* return prev because we just swapped (if not pruned) */
    free(probvec);
    return probvec_prev;
}
/* pruned_calc_prob_dist */


#ifdef PSEUDO_BINOMIAL
/* binomial test using poissbin. only good for high n and small prob.
 * returns -1 on error */
int
pseudo_binomial(long double *pvalue,
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
poissbin(long double *pvalue, const double *err_probs,
         const int num_err_probs, const int num_failures,
         const long long int bonf, const double sig)
{
    double *probvec = NULL;
    int errsv;
#if TIMING
    clock_t start = clock();
    int msec;
#endif
    *pvalue = LDBL_MAX;

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

    errno = 0;
    feclearexcept(FE_ALL_EXCEPT);

    *pvalue = expl(probvec[num_failures]); /* no need for tailsum here */

    errsv = errno;
    if (errsv || fetestexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW | FE_UNDERFLOW)) {
         if (*pvalue < DBL_EPSILON) {
              *pvalue = LDBL_MIN;/* to zero but prevent actual 0 value */
         } else {
              *pvalue = LDBL_MAX; /* otherwise set to 1 which might pass filters */
         }
    }

    return probvec;
}



/**
 * @brief
 *
 * pvalues computed for each of the NUM_NONCONS_BASES noncons_counts
 * will be written to snp_pvalues in the same order. If pvalue was not
 * computed (always insignificant) its value will be set to LDBL_MAX
 *
 */
int
snpcaller(long double *snp_pvalues,
          const double *err_probs, const int num_err_probs,
          const int *noncons_counts,
          const long long int bonf_factor, const double sig_level)
{
    double *probvec = NULL;
    int i;
    int max_noncons_count = 0;
    long double pvalue;

#if 0
    for (i=0; i<num_err_probs; i++) {
         fprintf(stderr,  "%f ", err_probs[i]);
    }
    fprintf(stderr,  "\n");
#endif

#ifdef DEBUG
    fprintf(stderr, "DEBUG(%s:%s():%d): num_err_probs=%d noncons_counts=%d,%d,%d bonf_factor=%lld sig_level=%f\n",
            __FILE__, __FUNCTION__, __LINE__,
            num_err_probs, noncons_counts[0], noncons_counts[1], noncons_counts[2],
            bonf_factor, sig_level);
#endif

    /* initialise empty results so that we can return anytime */
    for (i=0; i<NUM_NONCONS_BASES; i++) {
        snp_pvalues[i] = LDBL_MAX;
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

#if 0
    for (i=1; i<max_noncons_count+1; i++) {
        fprintf(stderr, "DEBUG(%s:%s():%d): prob for count %d=%Lg\n",
                __FILE__, __FUNCTION__, __LINE__,
                i, expl(probvec[i]));
    }
#endif

    if (pvalue * (double)bonf_factor > sig_level) {
#ifdef DEBUG
        fprintf(stderr, "DEBUG(%s:%s():%d): Most frequent SNV candidate already gets not signifcant pvalue of %Lg * %lld > %f\n",
                __FILE__, __FUNCTION__, __LINE__,
                pvalue, bonf_factor, sig_level);
#endif
        goto free_and_exit;
    }


    /* report p-value for each non-consensus base
     */
    for (i=0; i<NUM_NONCONS_BASES; i++) {
        if (0 != noncons_counts[i]) {
             int errsv;
             errno = 0;
             feclearexcept(FE_ALL_EXCEPT);

             pvalue = expl(probvec_tailsum(probvec, noncons_counts[i], max_noncons_count+1));

             errsv = errno;
             if (errsv || fetestexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW | FE_UNDERFLOW)) {
                  /* failed expl will set pvalue either to 0 or 1,
                   * both of which is not wanted here: this function
                   * should never return 0.0 but only just (LDBL_MIN)
                   * and 1.0 might vreate problems with high Bonf/Sig
                   * factors so we need to return a high value
                   * (LDBL_MAX)
                   */
                 if (pvalue < DBL_EPSILON) {
                       pvalue = LDBL_MIN;/* to zero but prevent actual 0 value */
                  } else {
                       pvalue = LDBL_MAX; /* otherwise set to 1 which might pass filters */
                  }
             }
            snp_pvalues[i] = pvalue;
#ifdef DEBUG
            fprintf(stderr, "DEBUG(%s:%s():%d): i=%d noncons_counts=%d max_noncons_count=%d pvalue=%Lg\n",
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
 * newer versions need the convoluted
 * gcc -Wall -g -std=gnu99 -O2 -DSNPCALLER_MAIN [-DUSE_SNPCALLER] -o snpcaller -I../uthash/ -I../libbam/ snpcaller.c utils.c log.c   plp.c samutils.c ../libbam/libbam.a -lm -lz -lpthread -DNDEBUG
 *

 Could use poibin for testing but parameter choice there is unclear

 library(poibin)
 # if pnorm is missing also do library(stats)

 pp=c(0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001)
 nerrs = 1
 # convert to success probabilities
 pp=1-pp
 > dpoibin(kk=length(pp)-nerrs, pp=pp)
 [1] 0.009910359
 > ppoibin(kk=length(pp)-nerrs, pp=pp)
 [1] 0.00995512
 # no approximation seems to work better:
 > ppoibin(kk=length(pp)-nerrs, pp=pp, method="NA")
 [1] 4.732391e-07

 ./snpcaller 10 1 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001num_trials=10 num_errs=1
 0.00896408

 ?!



 $ ./snpcaller 957 9  $(cat eprobs.txt)
 num_trials=957 num_errs=9
 1.77668e-05

 p = read.table('scratch/errprobs')
 pp = c(p)$V1
 nerrs = 9
 > ppoibin(kk=957-9, pp=1-pp)
 [1] 0.000262769
 > ppoibin(kk=957-9, pp=1-pp, method="NA")
 [1] 1.162356e-05

 ?!

 *
 */
int main(int argc, char *argv[]) {
     int num_trials;
     int num_errs;
     double *err_probs;
     int i;
     const float bonf = 1.0 ;
     const float sig = 1.0 ;

     verbose = 1;

     if (argc<4) {
          LOG_FATAL("%s\n", "need: num_trials num_errs p_e1 ... p_en");
          return -1;
     }

     num_trials = atoi(argv[1]);
     num_errs = atoi(argv[2]);
     if (argc-3 != num_trials) {
          LOG_FATAL("number of trials (%d) doesn't match number of error probabilities (%d)\n", num_trials, argc-3);
          exit(1);
     }
     err_probs = malloc(sizeof(double) * num_trials);
     for (i=3; i<argc; i++) {
          err_probs[i-3] = atof(argv[i]);
     }
     LOG_VERBOSE("num_trials=%d num_errs=%d\n", num_trials, num_errs);


#ifdef PSEUDO_BINOMIAL
     {
          if (-1 == pseudo_binomial(&pvalue,
                                    num_success, num_trials, succ_prob)) {
               LOG_ERROR("%s\n", "pseudo_binomial() failed");
               return -1;
          }
          printf("pseudo_binomial: %g\n", pvalue);
     }
#endif


#ifdef USE_SNPCALLER
     {
          long double snp_pvalues[NUM_NONCONS_BASES];
          int noncons_counts[NUM_NONCONS_BASES];
          noncons_counts[0] = num_errs;
          noncons_counts[1] = num_errs-1;
          noncons_counts[2] = num_errs-2;

          snpcaller(snp_pvalues, err_probs, num_trials, noncons_counts, bonf, sig);
          printf("prob from snpcaller(): (.. -2:%Lg .. -1:%Lg ..) = %Lg\n", snp_pvalues[2], snp_pvalues[1], snp_pvalues[0]);
     }
#else
     {
          double *probvec;
          long double pvalue;
          probvec = poissbin(&pvalue, err_probs, num_trials,
                             num_errs, bonf, sig);
          printf("%Lg\n", pvalue);
          free(probvec);
     }
#endif

     free(err_probs);
}
#endif
