/* -*- mode: c; tab-width: 4; c-basic-offset: 4;  indent-tabs-mode: nil -*- */

/*********************************************************************
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


#ifndef MAIN
/* must come first */
#include <Python.h>
#endif

#define TIMING 0

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include <ctype.h>
#include <float.h>

#include "cdflib.h"
#include "fet.h"

#if TIMING
#include <time.h>
#endif

#define PHREDQUAL_TO_PROB(phred) (pow(10.0, -1.0*(phred)/10.0))
#define PHREDQUAL_VALID_RANGE(phred) ((phred)>1 && (phred)<100)
#define PROB_TO_PHREDQUAL(prob) ((int)(-10.0 * log10(prob)))

#define MIN(X,Y) ((X) < (Y) ? (X) : (Y))
#define MAX(X,Y) ((X) > (Y) ? (X) : (Y))

#define LOGZERO -1e100 
#define FLOAT_DELTA 1e-32

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

#if 0
#define DEBUG_PY2C
#endif


int binom_sf(double *q, int num_trials, int num_successes,
             double prob_success);
double log_sum(double log_a, double log_b);
double log_diff(double log_a, double log_b);
double probvec_tailsum(double *probvec, int tail_startindex,
                       int probvec_len);
double *naive_calc_prob_dist(const int *quals, int N, int K);
double *pruned_calc_prob_dist(const int *quals, int N, int K, 
                      unsigned long int bonf_factor, double sig_level);
int snpcaller_qual(double *snp_pvalues, const int *phred_quals,
                   const int num_phred_quals, const int *noncons_counts,
                   const unsigned long int bonf_factor,
                   const double sig_level);
static PyObject *py_binom_sf(PyObject *self, PyObject *args);
static PyObject *py_binom_cdf(PyObject *self, PyObject *args);
static PyObject *py_kt_fisher_exact(PyObject *self, PyObject *args);
static PyObject *py_snpcaller_qual(PyObject *self, PyObject *args);
static PyObject *py_phredqual_to_prob(PyObject *self, PyObject *args);
static PyObject *py_prob_to_phredqual(PyObject *self, PyObject *args);





/**
 * @brief Survival function (1-cdf) computed by means of DCDFLIB to get rid of scipy
 * dependency. Same as scipy.stats.binom.sf().
 *
 * q will hold the probability of seeing num_successes or more in num_trials
 * given a probabilty of prob_success
 *
 * Returns non-zero status on failure
 *
 */
int binom_sf(double *q,
			 int num_trials, int num_successes, double prob_success) 
{
		int which=1;
		int status=1; /* error by default */
		double ompr = 1.0 - prob_success;
		double bound;
		double p;

		double s = (double)num_successes;
		double xn = (double)num_trials;
		double pr = (double)prob_success;

		(void) cdfbin(&which, &p, q,
			   &s, &xn, &pr, &ompr,
			   &status, &bound);

#if 0
		fprintf(stderr, "DEBUG(%s:%s:%d): in num_successes = %d\n", 
                __FILE__, __FUNCTION__, __LINE__, num_successes);
		fprintf(stderr, "DEBUG(%s:%s:%d): in num_trials = %d\n", 
                __FILE__, __FUNCTION__, __LINE__, num_trials);
		fprintf(stderr, "DEBUG(%s:%s:%d): in pr = %g\n", 
                __FILE__, __FUNCTION__, __LINE__, prob_success);
		fprintf(stderr, "DEBUG(%s:%s:%d): out p=%g\n", 
                __FILE__, __FUNCTION__, __LINE__, p);
		fprintf(stderr, "DEBUG(%s:%s:%d): out binom.sf() q=%g\n", 
                __FILE__, __FUNCTION__, __LINE__, *q);
		fprintf(stderr, "DEBUG(%s:%s:%d): out status=%d\n",
                __FILE__, __FUNCTION__, __LINE__, status);
		fprintf(stderr, "DEBUG(%s:%s:%d): out bound=%g\n", 
                __FILE__, __FUNCTION__, __LINE__, bound);
#endif
		
		return status;
}
/* end of binom_sf */



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
/* end of log_sum() */


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
/* end of log_diff() */



/**
 * @brief Computes sum of probvec values (log space) starting from (including)
 * tail_startindex to (excluding) probvec_len
 *
 */
double
probvec_tailsum(double *probvec, int tail_startindex, int probvec_len)
{
    double tailsum;
    int i;

    tailsum = probvec[tail_startindex];
    for (i=tail_startindex+1; i<probvec_len; i++) {
        tailsum = log_sum(tailsum, probvec[i]);
    }

    return tailsum;
}
/* end of probvec_tailsum() */


/**
 * @brief FIXME:missing-doc
 *
 *
 */
double *
naive_calc_prob_dist(const int *quals, int N, int K)
{
   double *probvec = NULL;
   double *probvec_prev = NULL;
   double *probvec_swp = NULL;

   int n;
   fprintf(stderr, "CRITICAL(%s:%s:%d): Use pruned_calc_prob_dist instead of me\n", 
           __FILE__, __FUNCTION__, __LINE__);


    if (NULL == (probvec = malloc((N+1) * sizeof(double)))) {
        fprintf(stderr, "FATAL: couldn't allocate memory at %s:%s():%d\n",
                __FILE__, __FUNCTION__, __LINE__);
        return NULL;
    }
    if (NULL == (probvec_prev = malloc((N+1) * sizeof(double)))) {
        fprintf(stderr, "FATAL: couldn't allocate memory at %s:%s():%d\n",
                __FILE__, __FUNCTION__, __LINE__);
        return NULL;
    }

    /* init */
    probvec_prev[0] = 0.0; /* 0.0 = log(1.0) */

    for (n=1; n<N+1; n++) {
        int k;
        double pn = PHREDQUAL_TO_PROB(quals[n-1]);
        double log_pn = log(pn);
        double log_1_pn = log1p(-pn);

        assert(PHREDQUAL_VALID_RANGE(quals[n-1]));

        k = 0;
        probvec[k] = probvec_prev[k] + log_1_pn;

        for (k=1; k<K; k++) {
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
/* end of naive_prob_dist */



/**
 * @brief FIXME:missing-doc
 *
 * Should really get rid of bonf_factor and sig_level here and
 * upstream as well
 *
 */
double *
pruned_calc_prob_dist(const int *quals, int N, int K, 
                      unsigned long int bonf_factor, double sig_level)
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
        double pn = PHREDQUAL_TO_PROB(quals[n-1]);
        double log_pn = log(pn);
        double log_1_pn = log1p(-pn); /* 0.0 = log(1.0) */
        
        /* test for valid phred quality boundaries */
        assert(PHREDQUAL_VALID_RANGE(quals[n-1]));

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
/* end of pruned_calc_prob_dist */




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
snpcaller_qual(double *snp_pvalues, 
               const int *phred_quals, const int num_phred_quals, 
               const int *noncons_counts, 
               const unsigned long int bonf_factor, const double sig_level)
{
    double *probvec = NULL;
    int i;
#if TIMING
    clock_t start = clock();
    int msec;
#endif
    int max_noncons_count = 0;
    double pvalue;

#ifdef DEBUG
    fprintf(stderr, "DEBUG(%s:%s():%d): num_phred_quals=%d noncons_counts=%d,%d,%d bonf_factor=%lu sig_level=%f\n", 
            __FILE__, __FUNCTION__, __LINE__, 
            num_phred_quals, noncons_counts[0], noncons_counts[1], noncons_counts[2],
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

#if TIMING
    start = clock();
#endif

#ifdef NAIVE
    probvec = naive_calc_prob_dist(phred_quals, num_phred_quals,
                                    max_noncons_count);
#else
    probvec = pruned_calc_prob_dist(phred_quals, num_phred_quals,
                                    max_noncons_count,
                                    bonf_factor, sig_level);
#endif

#if TIMING
    msec = (clock() - start) * 1000 / CLOCKS_PER_SEC;
    fprintf(stderr, "calc_prob_dist() took %d s %d ms\n", msec/1000, msec%1000);
#endif

    pvalue = exp(probvec[max_noncons_count]);
    if (pvalue * (double)bonf_factor >= sig_level) {
#ifdef DEBUG
        fprintf(stderr, "DEBUG(%s:%s():%d): Most frequent SNV candidate already gets not signifcant pvalue of %g * %ul >= %g\n", 
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
            pvalue = exp(probvec_tailsum(probvec, noncons_counts[i], max_noncons_count+1));
            snp_pvalues[i] = pvalue;
#ifdef DEBUG
            fprintf(stderr, "DEBUG(%s:%s():%d): pvalue=%g for noncons_counts[i]=%d bonf_factor=%ul\n", 
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
/* end of snpcaller_qual() */




#ifndef MAIN

/*
 * Python Doc
 * ----------
 *  
 * General:
 * - http://docs.python.org/extending/extending.html
 * - http://www.tutorialspoint.com/python/python_further_extensions.htm
 * 
 * Specific:
 * 
 * - Python C extension example:Translating a Python Sequence into a C Array with the PySequence_Fast Protocol:
 * http://book.opensourceproject.org.cn/lamp/python/pythoncook2/opensource/0596007973/pythoncook2-chp-17-sect-6.html
 * 
 * - Recipe 17.7. Accessing a Python Sequence Item-by-Item with the Iterator Protocol
 * http://book.opensourceproject.org.cn/lamp/python/pythoncook2/opensource/0596007973/pythoncook2-chp-17-sect-7.html
 * 
 *
 */



/**
 * @brief Python wrapper for binom_sf()
 *
 */
static PyObject *
py_binom_sf(PyObject *self, PyObject *args)
{
    int num_trials;
    int num_successes; 
    double prob_success;
    double out_prob;

    if (!PyArg_ParseTuple(args, "iid",
                          & num_successes,
                          & num_trials,
                          & prob_success)) {
        PyErr_SetString(PyExc_TypeError, "Failed to parse arguments.");
        return NULL;
    }

    if (num_successes<1) {
         return Py_BuildValue("d", 1.0);
    }

    /*if ((num_trials<num_successes) 
        || -> 1.0!
    if ((num_trials<0)
        ||
        (prob_success>1.0 || prob_success<0.0)
        ) {
        PyErr_SetString(PyExc_ValueError, "Arguments for binom_sf() don't make sense");
        return NULL;        
    }*/

    if (binom_sf(&out_prob, num_trials, num_successes, prob_success)) {
        PyErr_SetString(PyExc_TypeError, "binom_sf() failed");
        return NULL;
    }
    
    return Py_BuildValue("d", out_prob);
}
/* end of py_binom_sf() */


/**
 * @brief CDF 1.0-binom_sf()
 * 
 * @warn this is just a hack which might not have same precision as if
 * we used cdfbin directly
 * 
 */
static PyObject *
py_binom_cdf(PyObject *self, PyObject *args)
{
    int num_trials;
    int num_successes; 
    double prob_success;
    double out_prob;

    if (!PyArg_ParseTuple(args, "iid",
                          & num_successes,
                          & num_trials,
                          & prob_success)) {
        PyErr_SetString(PyExc_TypeError, "Failed to parse arguments.");
        return NULL;
    }

    /* mimick behaviour of scipy even though I don't understand their decisions */
    if (prob_success>1.0 || prob_success<0.0) {
        PyErr_SetString(PyExc_ValueError, "probability values passed to binom_cdf() don't make sense");
        return NULL;        
    }
    if (num_trials<0) {
         return Py_BuildValue("d", 1.0);
    }    
    if (num_successes<0) {
         return Py_BuildValue("d", 0.0);
    }
    if (num_successes>=num_trials) {
         return Py_BuildValue("d", 1.0);    
    }

    if (binom_sf(&out_prob, num_trials, num_successes, prob_success)) {
        PyErr_SetString(PyExc_TypeError, "binom_sf() failed");
        return NULL;
    }
    
    return Py_BuildValue("d", 1.0-out_prob);
}
/* end of py_binom_cdf() */



/**
 * @brief Python wrapper for kt_fisher_exact()
 *
 */
static PyObject *
py_kt_fisher_exact(PyObject *self, PyObject *args)
{
    double prob, left, right, twotail;
	int n11, n12, n21, n22;
    PyObject *py_pvalues;

    if (!PyArg_ParseTuple(args, "(ii)(ii)",
                          & n11, &n12, &n21, &n22)) {
        PyErr_SetString(PyExc_TypeError, "Failed to parse arguments.");
        return NULL;
    }

    if (n11<0 && n12<0 && n21<0 && n22<0) {
        PyErr_SetString(PyExc_ValueError, "Arguments cannot be negative");
        return NULL;
    }


    prob = kt_fisher_exact(n11, n12, n21, n22, &left, &right, &twotail);

    py_pvalues = Py_BuildValue("(ddd)", 
                               left, right, twotail);

    return py_pvalues;
}
/* end of py_kt_fisher_exact() */



/**
 * @brief Wrapper for snpcaller_qual()
 *
 */
static PyObject *
py_snpcaller_qual(PyObject *self, PyObject *args)
{
    /* in */
    PyObject *py_phred_quals;
    int *phred_quals;
    int num_phred_quals;
    int noncons_counts[NUM_NONCONS_BASES];
    unsigned long int bonf_factor;
    double sig_level;

    /* out */
    PyObject *py_snp_pvalues;
    double *snp_pvalues;

    /* aux */
    int i;

    if (NULL == (snp_pvalues = malloc(NUM_NONCONS_BASES * sizeof(double)))) {
        fprintf(stderr, "FATAL: couldn't allocate memory at %s:%s():%d\n",
                __FILE__, __FUNCTION__, __LINE__);
        return PyErr_NoMemory();    
    }
    

    if (!PyArg_ParseTuple(args, "O(iii)kd",
                          & py_phred_quals,
                          & noncons_counts[0],
                          & noncons_counts[1],
                          & noncons_counts[2],
                          & bonf_factor,
                          & sig_level)) {
        PyErr_SetString(PyExc_TypeError, "Failed to parse arguments.");
        return NULL;
    }
  

    /* parse py_phred_quals and save to phred_quals
     */
    py_phred_quals = PySequence_Fast(py_phred_quals,
                                     "argument must be iterable");
    if (! py_phred_quals) {
        return NULL;
    }
    num_phred_quals = PySequence_Fast_GET_SIZE(py_phred_quals);
    if (NULL == (phred_quals = malloc(num_phred_quals * sizeof(int)))) {
        Py_DECREF(py_phred_quals);
        fprintf(stderr, "FATAL: couldn't allocate memory at %s:%s():%d\n",
                __FILE__, __FUNCTION__, __LINE__);
        return PyErr_NoMemory();
    }
    for (i=0; i < num_phred_quals; i++) {
        PyObject *fitem;
        PyObject *item = PySequence_Fast_GET_ITEM(py_phred_quals, i);

        if (! item) {
            Py_DECREF(py_phred_quals);
            free(phred_quals);
            return NULL;
        }
        fitem = PyNumber_Int(item);
        if (! fitem) {
            Py_DECREF(py_phred_quals);
            free(phred_quals);
            PyErr_SetString(PyExc_TypeError, "All quality items must be numbers");
            return NULL;
        }
        phred_quals[i] = (int)PyInt_AsLong(fitem);
#ifdef DEBUG_PY2C
        fprintf(stderr, "DEBUG(%s:%s:%d): set phred_quals[%d]=%d\n",
                __FILE__, __FUNCTION__, __LINE__, i, phred_quals[i]);
#endif
        Py_DECREF(fitem);
    }
    Py_DECREF(py_phred_quals);



    if (snpcaller_qual(snp_pvalues, phred_quals, num_phred_quals, 
                  noncons_counts, bonf_factor, sig_level)) {
        PyErr_SetString(PyExc_TypeError, "snpcaller_qual() failed");
        py_snp_pvalues = NULL;
        goto free_and_return_py_snp_calues;
    }

    py_snp_pvalues = Py_BuildValue("(ddd)", 
                         snp_pvalues[0],
                         snp_pvalues[1],
                         snp_pvalues[2]);

 free_and_return_py_snp_calues:
    free(phred_quals);
    free(snp_pvalues);

    return py_snp_pvalues;
}
/* end of py_snpcaller_qual() */



/**
 * @brief Wrapper for PHREDQUAL_TO_PROB()
 *
 */
static PyObject *
py_phredqual_to_prob(PyObject *self, PyObject *args)
{
    int phredqual;

    if (!PyArg_ParseTuple(args, "i", &phredqual)) {
        PyErr_SetString(PyExc_TypeError, "Failed to parse arguments.");
        return NULL;
    }

#ifndef NDEBUG
    if (! (PHREDQUAL_VALID_RANGE(phredqual))) {
        PyErr_SetString(PyExc_TypeError, "phred quality out of valid boundaries");
        return NULL;
    }
#endif

    return Py_BuildValue("d", PHREDQUAL_TO_PROB(phredqual));
}
/* end of py_phredqual_to_prob() */



/**
 * @brief Wrapper for PROB_TO_PHREDQUAL
 *
 */
static PyObject *
py_prob_to_phredqual(PyObject *self, PyObject *args)
{
    int prob;

    fprintf(stderr, "CRITICAL(%s:%s:%d): untested function\n",
            __FILE__, __FUNCTION__, __LINE__);

    if (!PyArg_ParseTuple(args, "d", &prob)) {
        PyErr_SetString(PyExc_TypeError, "Failed to parse arguments.");
        return NULL;
    }

#ifndef NDEBUG
    if (prob<0.0 || prob>1.0) {
        PyErr_SetString(PyExc_TypeError, "prob out of valid boundaries");
        return NULL;
    }
#endif
    /* check for prob=0.0 */
    return Py_BuildValue("i", PROB_TO_PHREDQUAL(prob));
}
/* end of py_phredqual_to_prob() */


/**
 *
 */
static PyMethodDef
LoFreqExtMethods[] = {
    {"snpcaller_qual", py_snpcaller_qual, 
     METH_VARARGS, "Calculate SNPs for one pileup columns"},
    {"binom_sf", py_binom_sf, 
     METH_VARARGS, "Survival function (1-cdf)"},
    {"binom_cdf", py_binom_cdf, 
     METH_VARARGS, "Cumulative density function."},
    {"kt_fisher_exact", py_kt_fisher_exact,
     METH_VARARGS, "Fisher's Exact Test"},


    /* only for testing */
    {"_phredqual_to_prob", py_phredqual_to_prob,
     METH_VARARGS, "Convert Phred-quality to probability"},
    {"_prob_to_phredqual", py_prob_to_phredqual,
     METH_VARARGS, "Convert error-probability to Phred-quality"},

    /* sentinel */
    {NULL, NULL, 0, NULL}
};



/**
 * called on first import
 *
 */
PyMODINIT_FUNC
initlofreq_ext(void)
{
    (void) Py_InitModule("lofreq_ext", LoFreqExtMethods);
}
/* end of initlofreq_ext() */


#else


/**
 * @brief
 *
 * gcc -DMAIN -DDEBUG=1 -o lofreq_ext -ansi -Wall -pedantic  -g -Icdflib90/ cdflib90/dcdflib.c cdflib90/ipmpar.c lofreq_ext.c
 *
 */
int main(int argc, char **argv)
{
    double pvalues[3];
    unsigned long int bonf = 1;
    double sig = 0.05;

    if (0) {
        int num_phred_quals = 10;
        int phred_quals[] = {10, 10, 10, 10, 10, 10, 10, 10, 10, 10};
        int noncons_counts[3] = {0, 1, 2};
        snpcaller_qual(pvalues, 
                       phred_quals, (const int)num_phred_quals, (const int *)noncons_counts,
                       (const unsigned long int)bonf, (const double)sig);
    }

    if (0) {
        int num_phred_quals = 10;
        int phred_quals[] = {20, 20, 20, 20, 20, 20, 20, 20, 20, 20};
        int noncons_counts[3] = {0, 1, 4};
        snpcaller_qual(pvalues, 
                       phred_quals, (const int)num_phred_quals, (const int *)noncons_counts,
                       (const unsigned long int)bonf, (const double)sig);
    }

    if (1) {
        int num_phred_quals = 100;
        int phred_quals[] = {20, 20, 20, 20, 20, 20, 20, 20, 20, 20,
                             20, 20, 20, 20, 20, 20, 20, 20, 20, 20,
                             20, 20, 20, 20, 20, 20, 20, 20, 20, 20,
                             20, 20, 20, 20, 20, 20, 20, 20, 20, 20,
                             20, 20, 20, 20, 20, 20, 20, 20, 20, 20,
                             20, 20, 20, 20, 20, 20, 20, 20, 20, 20,
                             20, 20, 20, 20, 20, 20, 20, 20, 20, 20,
                             20, 20, 20, 20, 20, 20, 20, 20, 20, 20,
                             20, 20, 20, 20, 20, 20, 20, 20, 20, 20,
                             20, 20, 20, 20, 20, 20, 20, 20, 20, 20};
        int noncons_counts[3] = {1, 5, 10};
        snpcaller_qual(pvalues, 
                       phred_quals, (const int)num_phred_quals, (const int *)noncons_counts,
                       (const unsigned long int)bonf, (const double)sig);
    }


    if (0) {
        int num_phred_quals = 10;
        int phred_quals[] = {30, 30, 30, 30, 30, 30, 30, 30, 30, 30};
        int noncons_counts[3] = {0, 1, 4};
        snpcaller_qual(pvalues, 
                       phred_quals, (const int)num_phred_quals, (const int *)noncons_counts,
                       (const unsigned long int)bonf, (const double)sig);
    }

    if (0) {
        int phred_quals[10] = {30, 32, 32, 22, 25, 25, 25, 31, 25, 23};
        int num_phred_quals = 10;
        int noncons_counts[3] = {0, 1, 4};
        snpcaller_qual(pvalues, 
                       phred_quals, (const int)num_phred_quals, (const int *)noncons_counts,
                       (const unsigned long int)bonf, (const double)sig);
    }

    return EXIT_SUCCESS;
}
/* end of main() */


#endif
