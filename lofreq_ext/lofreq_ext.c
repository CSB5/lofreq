/* -*- mode: c; tab-width: 4; c-basic-offset: 4;  indent-tabs-mode: nil -*- */

#ifndef MAIN
/* must come first */
#include <Python.h>
#endif

#define TIMING 0
#define OPTIMIZE 1
#define DEFAULT_SIGLEVEL 0.05

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include <ctype.h>
#include "cdflib.h"

#if TIMING
#include <time.h>
#endif

/* Turn phred score into a probability */
#define PHREDQUAL_TO_PROB(pq) (pow(10.0, -1.0*(pq)/10.0))
#define PHREDQUAL_VALID_RANGE(pq) ((pq)>1 && (pq)<100)

/* Four nucleotides, with one consensus, makes three
       non-consensus */
#define NUM_NONCONS_BASES 3

#if 0
#define DEBUG
#endif

#if 0
#define DEBUG_PY2C
#endif






/**
 * @brief Survival function (1-cdf) computed by means of DCDFLIB to get rid of scipy
 * dependency. Same as scipy.stats.binom.sf().
 *
 * q will hold the probabilty of seeing num_successes or more in num_trials
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
                          & num_trials,
                          & num_successes,
                          & prob_success)) {
        PyErr_SetString(PyExc_TypeError, "Failed to parse arguments.");
        return NULL;
    }

    if (num_successes<1) {
         return Py_BuildValue("d", 1.0);
    }

    if ((num_trials<num_successes) 
        ||
        (num_trials<0) 
        ||
        (prob_success>1.0 || prob_success<0.0)
        ) {
        PyErr_SetString(PyExc_ValueError, "Arguments for binom_sf() don't make sense");
        return NULL;        
    }

    if (binom_sf(&out_prob, num_trials, num_successes, prob_success)) {
        PyErr_SetString(PyExc_TypeError, "binom_sf() failed");
        return NULL;
    }
    
    return Py_BuildValue("d", out_prob);
}
/* end of py_binom_sf() */




/**
 * @brief Computes log(exp(log_a) + exp(log_b))
 *
 * Taken from util.h of FAST source code:
 * http://www.cs.cornell.edu/~keich/FAST/fast.tar.gz
 *
 */
double
log_sum(double log_a, double log_b)
{
    if (log_a > log_b) {
        return log_a + log(1+exp(log_b-log_a));
    } else {
        return log_b + log(1+exp(log_a-log_b));
    }
}
/* end of log_sum() */



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
 * @brief Computes the following probability distribution for given for phred
 * scores:
 *       
 * P_{n}(X=k)  =   P_{n-1}(X=k) * (1 - p_{n})
 *                 + P_{n-1}(X=k-1) * p{_n}
 * 
 * Array quals has to be of total_num_bases size.
 * Non-consensus qualities should already be part of quals)
 * Values of returned vector are in log-space.
 * 
 * Caller has to free
 *
 */
double *
#if OPTIMIZE
calc_prob_dist(const int *quals, int total_num_bases, int max_noncons_count, 
               int bonf_factor, double sig_level)
#else
calc_prob_dist(const int *quals, int total_num_bases)
#endif
{
    double *probvec = NULL;
    double *probvec_prev = NULL;
    int n;


    if (NULL == (probvec = malloc((total_num_bases+1) * sizeof(double)))) {
        fprintf(stderr, "FATAL: couldn't allocate memory at %s:%s():%d\n",
                __FILE__, __FUNCTION__, __LINE__);
        return NULL;
    }
    if (NULL == (probvec_prev = malloc((total_num_bases+1) * sizeof(double)))) {
        fprintf(stderr, "FATAL: couldn't allocate memory at %s:%s():%d\n",
                __FILE__, __FUNCTION__, __LINE__);
        return NULL;
    }
#ifdef INIT_TO_FUNNY_VALUES
    for ... {
      probvec_prev[i] = 666.666;
      probvec[i] = 666.666;
    }
#endif


    /* init */
    probvec_prev[0] = log(1.0); /* = 0.0 */
    for (n=1; n<total_num_bases+1; n++) {
        int k;
        double pn = PHREDQUAL_TO_PROB(quals[n-1]);
        double log_pn = log(pn);
        double log_1_pn = log(1.0 - pn);
        /* test for valid phred quality boundaries */
        assert(PHREDQUAL_VALID_RANGE(quals[n-1]));

        k = 0;
        probvec[k] = probvec_prev[k] + log_1_pn;
        
        for (k=1; k<n; k++) {
            probvec[k] = log_sum(probvec_prev[k] + log_1_pn,
                                 probvec_prev[k-1] + log_pn);
        }
        
        k = n;
        probvec[k] = probvec_prev[k-1] + log_pn;

#if OPTIMIZE
        /* exit early if even the maximum number of noncons_counts are not
         * significant. a more elegant, faster way would be to
         * iterate from the end and keep the sum. but the modulo trick
         * here works as well and required less testing */
        if (n>max_noncons_count && n%100==0) {
            double pvalue = exp(probvec_tailsum(probvec, max_noncons_count, n));
            if (pvalue * (double)bonf_factor > sig_level) {
#if DEBUG
                fprintf(stderr,
                        "DEBUG(%s:%s:%d): early exit at n=%d with max_noncons_count=%d, bonf_factor=%d pvalue=%g sig_level=%f\n", 
                         __FILE__, __FUNCTION__, __LINE__,
                        n, max_noncons_count, pvalue,
                        bonf_factor, sig_level);
#endif
                probvec=NULL;
                break;
            }
        }
#endif

        memcpy(probvec_prev, probvec, (total_num_bases+1) * sizeof(double));
    }


#ifdef DEBUG
    for (n=1; n<total_num_bases+1; n++) {
        fprintf(stderr, "DEBUG(%s:%s:%d): probvec[%d]=%f\n", 
                __FILE__, __FUNCTION__, __LINE__, n, probvec[n]);
    }
#endif
    
    free(probvec_prev);
    
    return probvec;
}
/* end of calc_prob_dist() */





/**
 * @brief
 * 
 * driver function
 *
 * Call per pileup column
 *
 *
 * P-Values computed for each of the NUM_NONCONS_BASES noncons_counts
 * will be written to snp_pvalues in the same order.
 * P-Values are corrected by Bonferroni factor
 * 
 */
int
snpcaller_qual(double *snp_pvalues, 
               const int *phred_quals, const int num_phred_quals, 
               const int *noncons_counts, 
               const int bonf_factor, const double sig_level)
{
    double *probvec;
    int i;
#if TIMING
    clock_t start = clock();
    int msec;
#endif
    int max_noncons_count = 0;

#ifdef DEBUG
            fprintf(stderr, "DEBUG(%s:%s():%d): num_phred_quals=%d noncons_counts=%d bonf_factor=%d sig_level=%f\n", 
                    __FILE__, __FUNCTION__, __LINE__, 
                    num_phred_quals, noncons_counts, bonf_factor, sig_level);
#endif

    /* initialise empty results so that we can return anytime */
    for (i=0; i<NUM_NONCONS_BASES; i++) {
        snp_pvalues[i] = 1.0 * (double) bonf_factor;
    }

    for (i=0; i<NUM_NONCONS_BASES; i++) {
        if (noncons_counts[i] > max_noncons_count) {
            max_noncons_count = noncons_counts[i];
        }
    }

#if OPTIMIZE
    probvec = calc_prob_dist(phred_quals, num_phred_quals, max_noncons_count,
                             bonf_factor, sig_level);
#else
    probvec = calc_prob_dist(phred_quals, num_phred_quals);
#endif

#if TIMING
    msec = (clock() - start) * 1000 / CLOCKS_PER_SEC;
    fprintf(stderr, "calc_prob_dist() took %d s %d ms\n", msec/1000, msec%1000);
    start = clock();
#endif
    /* NULL means early exit */
    if (NULL == probvec) {
        return 0;
    }

    /* report p-value (sum of probabilities at end of distribution)
     * for each non-consensus base
     */
    for (i=0; i<NUM_NONCONS_BASES; i++) {        
        if (0 != noncons_counts[i]) {
            double tailsum;
            double pvalue;

            /* fast enough: no need to optimize */
            tailsum = probvec_tailsum(probvec, noncons_counts[i], num_phred_quals);
            pvalue = exp(tailsum);
            snp_pvalues[i] = pvalue * (double)bonf_factor;
#ifdef DEBUG
            fprintf(stderr, "DEBUG(%s:%s():%d): tailsum=%g pvalue=%g noncons_counts[i]=%d snp_pvalues[i]=%g bonf_factor=%d\n", 
                    __FILE__, __FUNCTION__, __LINE__, 
                    tailsum, pvalue, noncons_counts[i], snp_pvalues[i], bonf_factor);
#endif
        }
    }

    if (NULL != probvec) {
        free(probvec);
    }

    return 0;
}
/* end of snpcaller_qual() */




/* ---------------------------------------------------------------------- */
#if MAIN
/* ---------------------------------------------------------------------- */

/**
 * @brief
 *
 */
int main(int argc, char **argv)
{
#if 0
    char *bases = "CAAAATAAAN";
    int phred_quals[10] = {30, 32, 32, 18, 25, 25, 25, 31, 25, 23};
#endif

    char *bases = "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCACACCCCCCCCCCCCCCCCCCCCCACCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCTCCCCCCCCCCCCCCCTCCCCCCCCCCCCCCCCCCCCCCC";
    int phred_quals[152] = {35, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 35, 38, 37, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 37, 37, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 37, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 36, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 37, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 37, 38, 38, 37, 38, 38, 37, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 36};
    
    int cutoff = 20;

    snpcaller_qual(bases, phred_quals, cutoff, cutoff, DEFAULT_SIGLEVEL);

    return EXIT_SUCCESS;
}
/* end of main() */


/* ---------------------------------------------------------------------- */
#else
/* ---------------------------------------------------------------------- */


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
    int bonf_factor;
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
    

    if (!PyArg_ParseTuple(args, "O(iii)id",
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
 * @brief Wrapper for probvec_tailsum()
 * 
 */
static PyObject *
py_probvec_tailsum(PyObject *self, PyObject *args)
{
    PyObject *py_probvec = NULL; /* input */
    double *probvec = NULL;
    int probvec_len = 0;
    int tail_startindex = 0;
    double tailsum = 0.0; /* value to compute/return */
    int i = 0;


	if (! PyArg_ParseTuple(args, "Oi", &py_probvec, &tail_startindex)) {
        PyErr_SetString(PyExc_TypeError, "Failed to parse arguments.");
        return NULL;
    }

    py_probvec = PySequence_Fast(py_probvec, "argument must be iterable");
    if (! py_probvec) {
        return NULL;
    }

    /* parse list py_probvec (of size probvec_len) to array of doubles
     * probvec
     */
    probvec_len = PySequence_Fast_GET_SIZE(py_probvec);
    if (NULL == (probvec = malloc(probvec_len * sizeof(double)))) {
        Py_DECREF(py_probvec);
        
        fprintf(stderr, "FATAL: couldn't allocate memory at %s:%s():%d\n",
                __FILE__, __FUNCTION__, __LINE__);
        return PyErr_NoMemory();
    }
    for (i=0; i < probvec_len; i++) {
        PyObject *fitem;
        PyObject *item = PySequence_Fast_GET_ITEM(py_probvec, i);

        if (!item) {
            Py_DECREF(py_probvec);
            free(probvec);
            return NULL;
        }
        fitem = PyNumber_Float(item);
        if (!fitem) {
            Py_DECREF(py_probvec);
            free(probvec);
            PyErr_SetString(PyExc_TypeError, "all items must be numbers");
            return NULL;
        }
        probvec[i] = PyFloat_AS_DOUBLE(fitem);
#ifdef DEBUG_PY2C
        fprintf(stderr, "DEBUG(%s:%s:%d): set probvec[%d]=%f\n",
                __FILE__, __FUNCTION__, __LINE__, i, probvec[i]);
#endif
        Py_DECREF(fitem);
    }
    Py_DECREF(py_probvec);

    
   /* compute tailsum and return
    */
    tailsum = probvec_tailsum(probvec, tail_startindex, probvec_len);
    free(probvec);
    return Py_BuildValue("d", tailsum);
}
/* end of py_probvec_tailsum() */




/**
 * @brief Wrapper for calc_prob_dist()
 *
 */
static PyObject *
py_calc_prob_dist(PyObject *self, PyObject *args)
{
    PyObject* py_phred_quals = NULL; /* handed down as list */
    PyObject *py_probvec = NULL; /* returned */
    int num_quals = 0;
    int *phred_quals = NULL;
    double *probvec = NULL;
    int i = 0;


	if (! PyArg_ParseTuple(args, "O", &py_phred_quals)) {
        PyErr_SetString(PyExc_TypeError, "Failed to parse arguments.");
        return NULL;
    }

    py_phred_quals = PySequence_Fast(py_phred_quals, "argument must be iterable");
    if (! py_phred_quals) {
        return NULL;
    }


    /* parse list py_phred_quals (of size num_quals) to array of ints
     * phred_quals
     *
     */
    num_quals = PySequence_Fast_GET_SIZE(py_phred_quals);
    if (NULL == (phred_quals = malloc(num_quals * sizeof(int)))) {
        Py_DECREF(py_phred_quals);
        
        fprintf(stderr, "FATAL: couldn't allocate memory at %s:%s():%d\n",
                __FILE__, __FUNCTION__, __LINE__);
        return PyErr_NoMemory( );
    }
    for (i=0; i < num_quals; i++) {
        PyObject *fitem;
        PyObject *item = PySequence_Fast_GET_ITEM(py_phred_quals, i);

        if (!item) {
            Py_DECREF(py_phred_quals);
            free(phred_quals);
            return NULL;
        }
        fitem = PyNumber_Int(item);
        if (!fitem) {
            Py_DECREF(py_phred_quals);
            free(phred_quals);
            PyErr_SetString(PyExc_TypeError, "all items must be numbers");
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


    /* compute probabilities and build result 
     */
#if OPTIMIZE
    fprintf(stderr, "WARNING: %s not changed to new API yet\n", __FUNCTION__);
    exit(1);
#else
    probvec = calc_prob_dist(phred_quals, num_quals);
#endif

    if (! probvec) {
        free(phred_quals);
        return NULL;
    }

    /* FIXME: this apparently leaks according to
     * http://stackoverflow.com/questions/3512414/does-this-pylist-appendlist-py-buildvalue-leak
     * so use PyList_New(num_quals+1) and PyList_SET_ITEM 
     *
     */
    py_probvec = PyList_New(0);
    for (i=0; i<=num_quals; i++) {
        PyList_Append(py_probvec, PyFloat_FromDouble(probvec[i]));
    }

    free(phred_quals);
    free(probvec);

    return py_probvec;
}
/* end of py_calc_prob_dist() */



/**
 *
 */
static PyMethodDef
LoFreqExtMethods[] = {
    /* main function */
    {"snpcaller_qual", py_snpcaller_qual, 
     METH_VARARGS, "Calculate SNPs for one pileup columns"},

    {"binom_sf", py_binom_sf, 
     METH_VARARGS, "Survival function (1-cdf)"},

    /* only for testing */
    {"_calc_prob_dist", py_calc_prob_dist, 
     METH_VARARGS, "Compute probability distribution"},
    {"_phredqual_to_prob", py_phredqual_to_prob,
     METH_VARARGS, "Convert PHRED quality to probability"},
    {"_probvec_tailsum", py_probvec_tailsum,
     METH_VARARGS, "Computes sum of probability distribtuion tail"},
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


#endif
