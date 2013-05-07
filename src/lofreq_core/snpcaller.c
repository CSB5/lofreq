/* -*- c-file-style: "k&r"; indent-tabs-mode: nil; -*-
 *
 *
 * Copyright (C) 2011, 2012 Genome Institute of Singapore
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the GNU
 * General Public License for more details.
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


#ifndef MIN
#define MIN(X,Y) ((X) < (Y) ? (X) : (Y))
#endif
#ifndef MAX
#define MAX(X,Y) ((X) > (Y) ? (X) : (Y))
#endif

#define LOGZERO -1e100 
/* FIXME shouldn't we use something from float.h ? */

/* Four nucleotides, with one consensus, makes three
   non-consensus bases */
#define NUM_NONCONS_BASES 3

#if 0
#define DEBUG
#endif

#if 0
#define TRACE
#endif

#if 0
#define NAIVE
#endif


double log_sum(double log_a, double log_b);
double log_diff(double log_a, double log_b);
double probvec_tailsum(const double *probvec, int tail_startindex,
                       int probvec_len);
double *naive_calc_prob_dist(const double *err_probs, int N, int K);
double *pruned_calc_prob_dist(const double *err_probs, int N, int K, 
                      long long int bonf_factor, double sig_level);


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
 * @brief FIXME:missing-doc
 *
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
             /* FIXME clang: The left operand of '+' is a garbage value */
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
 * @brief FIXME:missing-doc
 *
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
        return NULL;
    }

    for (n=0; n<N; n++) {
         assert(err_probs[n]>=0.0 && err_probs[n]<=1.0);
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

        assert(! isinf(probvec[0])); /* used to happen when first q=0 */

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
#ifdef TRACE
		fprintf(stderr, "DEBUG(%s:%s:%d): setting probvec_prev[n=%d] to LOGZERO\n", 
                __FILE__, __FUNCTION__, __LINE__, n);
#endif
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
        for (k=0; k<=MIN(n,K-1); k++) {
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
            assert(probvec_prev[K]<=0.0 && probvec_prev[K-1]<=0.0);
            probvec[K] = log_sum(probvec_prev[K], probvec_prev[K-1]+log_pn);
            pvalue = exp(probvec[K]);
#ifdef TRACE
            fprintf(stderr, "DEBUG(%s:%s:%d): pvalue=%g (probvec[K=%d]=%g) at n=%d\n", 
                    __FILE__, __FUNCTION__, __LINE__, pvalue, K, probvec[K], n);
#endif

            if (pvalue * (double)bonf_factor >= sig_level) {
#ifdef DEBUG
		fprintf(stderr, "DEBUG(%s:%s:%d): early exit at n=%d with pvalue %g\n", 
                __FILE__, __FUNCTION__, __LINE__, n, pvalue);
#endif
               goto free_and_exit;
            }
        }

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


#ifdef RETURNS_WRONG_PVALUES_FOR_HIGH_PROBS
/* binomial test using poissbin. inefficient but allows us to reuse
 * existing code. returns -1 on error */
int
pseudo_binomial(double *pvalue, 
                int num_success, int num_trials, double succ_prob)
{
     const long long int bonf = 1.0;
     const double sig = 1.0;
     /* returned pvalues > sig/bonf are not computed properly */
     double *probvec = NULL;
     double *probs;
     int i;

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
 *  to check pvalues for other numbers < (original num_failures), like
 *  so: exp(probvec_tailsum(probvec, smaller_numl, orig_num+1)) but
 *  only if first pvalue was below limits implied by bonf and sig.
 *  default pvalue is DBL_MAX (1 might still be significant).
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

    *pvalue = exp(probvec[num_failures]);
    assert(! isnan(*pvalue));
    return probvec;
}



/**
 * @brief
 * 
 * Call per pileup column
 *
 * P-Values computed for each of the NUM_NONCONS_BASES noncons_counts
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
#ifdef TRACE
    for (i=1; i<max_noncons_count+1; i++) {        
        fprintf(stderr, "DEBUG(%s:%s():%d): pvalue=%g for noncons_counts %d\n", 
                __FILE__, __FUNCTION__, __LINE__, 
                exp(probvec_tailsum(probvec, i, max_noncons_count+1)), i);
    }
#endif
    for (i=0; i<NUM_NONCONS_BASES; i++) { 
        if (0 != noncons_counts[i]) {
#ifdef DEBUG
                  fprintf(stderr, "DEBUG(%s:%s():%d): i=%d noncons_counts=%d max_noncons_count=%d probvec_tailsum=%g\n", 
                          __FILE__, __FUNCTION__, __LINE__, 
                          i, noncons_counts[i], max_noncons_count, probvec_tailsum(probvec, noncons_counts[i], max_noncons_count+1));                  
#endif
            pvalue = exp(probvec_tailsum(probvec, noncons_counts[i], max_noncons_count+1));
            snp_pvalues[i] = pvalue;
#ifdef DEBUG
            fprintf(stderr, "DEBUG(%s:%s():%d): pvalue=%g for noncons_counts[i]=%d bonf_factor=%lld\n", 
                    __FILE__, __FUNCTION__, __LINE__, 
                    pvalue, noncons_counts[i], bonf_factor);
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
gcc -pedantic -Wall -g -std=gnu99 -O2 -DSNPCALLER_MAIN -o snpcaller snpcaller.c utils.c log.c
*/

int main(int argc, char *argv[]) {
     int num_success;
     int num_trials;
     double succ_prob;
     double pvalue;
     verbose = 1;
     if (argc<4) {
          LOG_ERROR("%s\n", "need num_success num_trials and succ_prob as args");
          return -1;
     }

     num_success = atoi(argv[1]);
     num_trials = atoi(argv[2]);
     succ_prob = atof(argv[3]);


     LOG_VERBOSE("num_success=%d num_trials=%d succ_prob=%f\n", num_success, num_trials, succ_prob);
     if (-1 == pseudo_binomial(&pvalue, 
                               num_success, num_trials, succ_prob)) {
          LOG_ERROR("%s\n", "pseudo_binomial() failed");
          return -1;
     }

     printf("%g\n", pvalue);
}
#endif
