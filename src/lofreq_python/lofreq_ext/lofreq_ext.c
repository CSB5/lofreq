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

/* Python API for lofreq_core C functions
 *
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

/* must come first: */
#include <Python.h> 

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
#include "binom.h"
#include "snpcaller.h"
#include "bam2depth.h"
#include "sam_view.h"
#include "utils.h"

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


static PyObject *py_binom_sf(PyObject *self, PyObject *args);
static PyObject *py_binom_cdf(PyObject *self, PyObject *args);
static PyObject *py_kt_fisher_exact(PyObject *self, PyObject *args);
static PyObject *py_snpcaller(PyObject *self, PyObject *args);
static PyObject *py_phredqual_to_prob(PyObject *self, PyObject *args);
static PyObject *py_prob_to_phredqual(PyObject *self, PyObject *args);



/**
 * @brief Python wrapper for depth_stats
 *
 */
static PyObject *
py_depth_stats(PyObject *self, PyObject *args, PyObject *keywds) {
    double depth_mean;
    long int depth_num_nonzero_pos;
    static char *kwlist[] = {
        "bam", "region", "bed", "min_baseq", "min_mapq", NULL};
    char *bam = NULL;
    char *region = NULL;
    char *bed = NULL;
    int min_baseq = 0;
    int min_mapq = 0;
    
    if (!PyArg_ParseTupleAndKeywords(
            args, keywds, "s|zzii", kwlist, 
            &bam, &region, &bed, &min_baseq, &min_mapq)) {
        PyErr_SetString(PyExc_TypeError, "Failed to parse arguments.");
        return NULL; 
    }

#if 0
    fprintf(stderr, "DEBUG bam=%s region=%s bed=%s min_baseq=%d min_mapq=%d",
           bam, region, bed, min_baseq, min_mapq);
#endif

    if (0 != depth_stats(&depth_mean, &depth_num_nonzero_pos,
                         bam, region, bed, &min_baseq, &min_mapq)) {
        PyErr_SetString(PyExc_TypeError, "depth_stats failed");
        return NULL;
    }

    return Py_BuildValue("(dl)", depth_mean, depth_num_nonzero_pos);
}
/* end of py_depth_stats() */



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
 * @brief Python wrapper for getting header from a SAM/BAM file
 *
 */
static PyObject *
py_read_sam_header(PyObject *self, PyObject *args)
{
    int is_bam = 1;
    char *bam_file = NULL;
    sam_view_opts_t *svo = NULL;
	char *header_buf = NULL; 
    int header_size = 0;
    PyObject *py_header;

    if (!PyArg_ParseTuple(args, "s|i",
                          & bam_file,
                          & is_bam)) {
        PyErr_SetString(PyExc_TypeError, "Failed to parse arguments.");
        return NULL;
    }

    if (new_sam_view_opts(&svo)) {
        PyErr_SetString(PyExc_TypeError, "new_sam_view_opts() failed\n");
        return NULL;
    }
    if (! is_bam) {
        svo->is_bamin = 0;
    }
    svo->fn_out = strdup(tmpnam(NULL));
    svo->is_header_only = 1;
    
    if (main_samview(bam_file, svo)) {
        PyErr_SetString(PyExc_TypeError, "main_samview() failed\n");
        free_sam_view_opts(&svo);
        return NULL;
    }

    header_size = ae_load_file_to_memory(svo->fn_out, &header_buf);
    if (header_size < 0) {
        PyErr_SetString(PyExc_TypeError, "Couldn't read header file\n");
        free_sam_view_opts(&svo);
        return NULL;        
    }

    py_header = Py_BuildValue("s", header_buf);

    unlink(svo->fn_out);
    free_sam_view_opts(&svo);
    free(header_buf);

    return py_header;
}
/* end of py_read_sam_header() */


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
 * @brief Wrapper for snpcaller()
 *
 */
static PyObject *
py_snpcaller(PyObject *self, PyObject *args)
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



    if (snpcaller(snp_pvalues, phred_quals, num_phred_quals, 
                  noncons_counts, bonf_factor, sig_level)) {
        PyErr_SetString(PyExc_TypeError, "snpcaller() failed");
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
/* end of py_snpcaller() */



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
    {"snpcaller_qual", py_snpcaller, 
     METH_VARARGS, "Calculate SNPs for one pileup columns"},
    {"binom_sf", py_binom_sf, 
     METH_VARARGS, "Survival function (1-cdf)"},
    {"binom_cdf", py_binom_cdf, 
     METH_VARARGS, "Cumulative density function."},
    {"kt_fisher_exact", py_kt_fisher_exact,
     METH_VARARGS, "Fisher's Exact Test"},
    {"read_sam_header", py_read_sam_header,
     METH_VARARGS, "Get header from SAM/BAM file"},
    {"depth_stats", (PyCFunction)py_depth_stats,
     METH_VARARGS|METH_KEYWORDS, "BAM depth stats"},

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


