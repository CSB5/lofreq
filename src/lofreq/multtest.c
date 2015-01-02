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
#include <math.h>
#include <string.h>
#include <ctype.h>

/* lofreq includes */
#include "utils.h"
#include "multtest.h"
#ifdef MULTTEST_TEST
#include "time.h"
#endif


typedef struct {
     double p;
     long int i;
} ixp_t; /* indexed p-value */

int
ixp_dbl_cmp(const void *a, const void *b)
{
     const ixp_t *ia = (const ixp_t *)a;
     const ixp_t *ib = (const ixp_t *)b;
     return dbl_cmp(&ia->p, &ib->p);
}


/* writes corrected values to given data array 
 * 
 * data is of size size. alpha = type-1 error cutoff for each test
 *
 * will use size as bonferroni factor if num_test (AKA bonf fac) < 1
 */
void
bonf_corr(double data[], long int size, long int num_tests)
{
     long int i;
     long int bonf_fac;

     if (num_tests<1) {
          bonf_fac = size;
     } else {
          bonf_fac = num_tests;
     }

     for (i=0; i<size; i++) {
          data[i] *= bonf_fac;
     }
}


/* writes corrected values to array
 *
 * data is of size size. alpha = type-1 error cutoff for each test 
 *
 * NOTE: only values that were originally below alpha are corrected!
 * 
 */
void
holm_bonf_corr(double data[], long int size, double alpha, long int num_tests)
{
     long int i;
     long int lp = size;
     double tp;
     double pp;
     ixp_t *iarr;

     iarr = malloc(size * sizeof(ixp_t));
     
     if (num_tests<1) {
          lp = size;
     } else {
          lp = num_tests;
     }

     /* first index the pvalues and store them in the ixp data */
     for (i = 0; i < size; i++) {
          iarr[i].i = i;
          iarr[i].p = data[i];
     }
     qsort(iarr, size, sizeof(ixp_t), ixp_dbl_cmp);

     pp = iarr[0].p;
     for (i = 0; i < size; i++) {
          /* if the pvalue is different update lp according to how
           * many pvalues have been seen and update the previous
           * pvalue seen */
          if (dbl_cmp(&iarr[i].p, &pp) != 0) {
               /* lp = size - i; */
               if (num_tests<1) {
                    lp = size-i;
               } else {
                    lp = num_tests-i;
               }

               pp = iarr[i].p;
          }
          /* if below alpha, correct the pvalue in the original
           * data */
          tp = iarr[i].p * 1. / lp;
          if (dbl_cmp(&tp, &alpha) < 0) {
               data[iarr[i].i] = iarr[i].p * lp;
          }
     }
     free(iarr);
}


/* will use size instead of num_tests if num_tests>0.
 *
 * irejected is a pointer to a 1D-array of indices of rejected (i..e
 * significant values). It is allocated here, i.e. user has to free.
 *
 * content of data will not be overwritten
 */
long int
fdr(double data[], long int size, double alpha, long int num_tests, long int **irejected)
{

     ixp_t *iarr;
     long int i;
     long int nrejected = 0;
     long int n;

     iarr = malloc(size * sizeof(ixp_t));
     if (num_tests<1) {
          n = size;
     } else {
          n = num_tests;
     }
     /* first index the pvalues and store them in ixp data to sort indices*/
     for (i = 0; i < size; i++) {
          iarr[i].i = i;
          iarr[i].p = data[i];
     }
     qsort(iarr, size, sizeof(ixp_t), ixp_dbl_cmp);

     /* starting from the largest rank, evaluate p(m) < alpha * (m/M) where m is the
      * rank and M is the total number of pvalues. If true, reject pvals 1..m */
     for (i = size; i > 0; i--) { /* ranks are 1-based */
          if (iarr[i-1].p < (alpha*i/(float)n)) {
               nrejected = i; /* therefore, nrejected includes first rejected (0-based) */
               break;
          }
     }

     /* return data of indices to rejected pvalues */
     *irejected = NULL;
     if (nrejected) {
          (*irejected) = (long int*) malloc(nrejected * sizeof(long int));
          for (i = 0; i < nrejected; i++) {
               /* printf("%d\t%f\t%d\n", iarr[i].i, iarr[i].p); */
               (*irejected)[i] = iarr[i].i;
          }
     }
     free(iarr);
     return nrejected;
}

int
mtc_str_to_type(char *t) {
     if (0 == strcmp(t, "bonf") || 0 == strcmp(t, "bonferroni")) {
          return MTC_BONF;
     } else if (0 == strcmp(t, "holm") || 0 == strcmp(t, "holmbonf") ||  0 == strcmp(t, "holm-bonf") || 0 == strcmp(t, "holmbonferroni")) {
          return MTC_HOLMBONF;
     } else if (0 == strcmp(t, "fdr")) {
          return MTC_FDR;
     } else {
          return -1;
     }
}


void
mtc_str(char *buf, int mtc_type) {
     char *str;
     int i;
     str = mtc_type_str[mtc_type];
     strcpy(buf, &(str[4]));
     for (i=0; i<strlen(buf); i++) {
          buf[i] = tolower(buf[i]);
     }
}


/* stand alone test

# R results as reference:
> p = c(2.354054e-07,2.101590e-05,2.576842e-05,9.814783e-05,1.052610e-04,1.241481e-04,1.325988e-04,1.568503e-04,2.254557e-04,3.795380e-04,6.114943e-04,1.613954e-03,3.302430e-03,3.538342e-03,5.236997e-03,6.831909e-03,7.059226e-03,8.805129e-03,9.401040e-03,1.129798e-02,2.115017e-02,4.922736e-02,6.053298e-02,6.262239e-02,7.395153e-02,8.281103e-02,8.633331e-02,1.190654e-01,1.890796e-01,2.058494e-01,2.209214e-01,2.856000e-01,3.048895e-01,4.660682e-01,4.830809e-01,4.921755e-01,5.319453e-01,5.751550e-01,5.783195e-01,6.185894e-01,6.363620e-01,6.448587e-01,6.558414e-01,6.885884e-01,7.189864e-01,8.179539e-01,8.274487e-01,8.971300e-01,9.118680e-01,9.437890e-01)
> sum(p < 0.05)
[1] 22
> sum(p.adjust(p, "BH") < 0.05)
[1] 20
> sum(p.adjust(p, "BH", 1000) < 0.05)
[1] 10
> sum(p.adjust(p, "BH", 100) < 0.001)
[1] 3
sum(p.adjust(p, "BH", 10000) < 1)
[1] 11

# Results from standalone binary:
ps="2.354054e-07 2.101590e-05 2.576842e-05 9.814783e-05 1.052610e-04  1.241481e-04 1.325988e-04 1.568503e-04 2.254557e-04 3.795380e-04 6.114943e-04 1.613954e-03 3.302430e-03 3.538342e-03 5.236997e-03 6.831909e-03 7.059226e-03 8.805129e-03 9.401040e-03 1.129798e-02 2.115017e-02 4.922736e-02 6.053298e-02 6.262239e-02 7.395153e-02 8.281103e-02 8.633331e-02 1.190654e-01 1.890796e-01 2.058494e-01 2.209214e-01 2.856000e-01 3.048895e-01 4.660682e-01 4.830809e-01 4.921755e-01 5.319453e-01 5.751550e-01 5.783195e-01 6.185894e-01 6.363620e-01 6.448587e-01 6.558414e-01 6.885884e-01 7.189864e-01 8.179539e-01 8.274487e-01 8.971300e-01 9.118680e-01 9.437890e-01"
$./multtest 50 0.05 $ps
 20 rejected with alpha 0.050000 and 50 tests
$ ./multtest 1000 0.05 $ps
 10 rejected with alpha 0.050000 and 1000 tests
$ ./multtest 100 0.001 $ps
 3 rejected with alpha 0.001000 and 100 tests
$ ./multtest 10000 1 $ps
 11 rejected with alpha 1.000000 and 10000 tests


 gcc -o multtest multtest.c utils.c log.c -ansi -Wall -DMULTTEST_STANDALONE -I../uthash/
*/
#ifdef MULTTEST_STANDALONE
int main(int argc, char *argv[])
{
     int i;
     int ntests;
     float alpha;
     double *data;
     int data_size;
     int nrejected;
     long int *irejected;/* indices of rejected i.e. significant values */


     if (argc<4) {
          fprintf(stderr, "Usage: %s numtests alpha p1 ... pn\n", argv[0]);
          exit(1);
     }
     /*fprintf(stderr, "argc=%d\n", argc);*/
     ntests = atoi(argv[1]);
     alpha = atof(argv[2]);
     fprintf(stderr, "ntests=%d alpha=%f\n", ntests, alpha);
     data_size = argc-3;
     if (ntests<data_size) {
          fprintf(stderr, "FATAL: ntests=%d < data_size=%d\n", ntests, data_size);
          exit(1);
     }
     data = malloc((data_size) * sizeof(double));
     for (i=0; i<data_size; i++) {
          data[i] = atof(argv[i+3]);
          fprintf(stderr, "DEBUG data[%d]=%f\n", i, data[i]);
     }

     nrejected = fdr(data, data_size, alpha, ntests, &irejected);
     
     printf ("%d rejected with alpha %f and %d tests: ", nrejected, alpha, ntests);
     if (nrejected) {
          for (i = 0; i < nrejected; i++) {
               printf("%f, ", data[irejected[i]]);
          }
     } else {
          printf("None");
     }
     printf ("\n\n");
     
     free(data);
     free(irejected);
     exit(0);
}
#endif


/* gcc -o multtest multtest.c utils.c log.c -ansi -Wall -DMULTTEST_TEST */
#ifdef MULTTEST_TEST


/* http://stackoverflow.com/questions/6127503/shuffle-array-in-c
 * answer by John Leehey
 *
 * arrange the N elements of ARRAY in random order.
 * Only effective if N is much smaller than RAND_MAX;
 * if this may not be the case, use a better random
 * number generator. */
#define NELEMS(x)  (sizeof(x) / sizeof(x[0]))
static void shuffle(void *array, size_t n, size_t size) {
    char tmp[size];
    char *arr = array;
    size_t stride = size * sizeof(char);

    if (n > 1) {
        size_t i;
        for (i = 0; i < n - 1; ++i) {
            size_t rnd = (size_t) rand();
            size_t j = i + rnd / (RAND_MAX / (n - i) + 1);

            memcpy(tmp, arr + j * stride, size);
            memcpy(arr + j * stride, arr + i * stride, size);
            memcpy(arr + i * stride, tmp, size);
        }
    }
}

int main(int argc, char *argv[])
{
     int i;
     int ntests;
     srand(time(NULL));


     {
          double data[] = {0.000001, 0.008, 0.039, 0.041, 0.042, 0.06, 0.074, 0.205, 0.212, 0.216, 0.222, 0.251, 0.269, 0.275, 0.34, 0.341, 0.384, 0.569, 0.594, 0.696, 0.762, 0.94, 0.942, 0.975, 0.986};
          int data_len = 25;
          float alpha = 0.25;
          int num_tests = data_len;
          long int* irejected;/* indices of rejected i.e. significant values */
          int nrejected;
          int exp_sig = 5;

          /* FIXME use also for bonferroni */
          printf("*** FDR test with data from http://www.biostathandbook.com/multiplecomparisons.html\n\n");

          nrejected = fdr(data, data_len, alpha, num_tests, &irejected);
          printf ("original:\t");
          if (nrejected) {
               for (i = 0; i < nrejected; i++) {
                    printf("%f, ", data[irejected[i]]);
               }
          } else {
               printf("None");
          }
          printf ("\n");
                     
          if (nrejected==exp_sig) {/* FIXME should also test values */
               printf("PASS\n\n");
          } else {
               printf("FAIL\n\n");
               exit(1);
          }
          free(irejected);



          int cap = 10;
          printf ("capped to %d:\t", cap);
          nrejected = fdr(data, cap, alpha, num_tests, &irejected);
          if (nrejected) {
               for (i = 0; i < nrejected; i++) {
                    printf("%f, ", data[irejected[i]]);
               }
          } else {
               printf("None");
          }
          printf ("\n");
                     
          if (nrejected==exp_sig) {/* FIXME should also test values */
               printf("PASS\n\n");
          } else {
               printf("FAIL\n\n");
               exit(1);
          }
          free(irejected);


          printf ("shuffled:\t");
          shuffle(data, NELEMS(data), sizeof(data[0]));
          nrejected = fdr(data, data_len, alpha, num_tests, &irejected);
          if (nrejected) {
               for (i = 0; i < nrejected; i++) {
                    printf("%f, ", data[irejected[i]]);
               }
          } else {
               printf("None");
          }
          printf ("\n");
                     
          if (nrejected==exp_sig) {/* FIXME should also test values */
               printf("PASS\n\n");
          } else {
               printf("FAIL\n\n");
               exit(1);
          }
          free(irejected);
     }

     {
          double data[] = {0.010,  0.013, 0.014, 0.190, 0.350, 0.500, 0.630, 0.670, 0.750, 0.810};
          int data_len = 10;
          float alpha = 0.05;
          int num_tests = data_len;
          long int* irejected;/* indices of rejected i.e. significant values */
          int nrejected;
          int exp_sig = 3;

          printf("*** FDR test with data from http://www.unc.edu/courses/2007spring/biol/145/001/docs/lectures/Nov12.html\n\n");

          nrejected = fdr(data, data_len, alpha, num_tests, &irejected);
          printf ("original:\t");
          if (nrejected) {
               for (i = 0; i < nrejected; i++) {
                    printf("%f, ", data[irejected[i]]);
               }
          } else {
               printf("None");
          }
          printf ("\n");
                     
          if (nrejected==exp_sig) {/* FIXME should also test values */
               printf("PASS\n\n");
          } else {
               printf("FAIL\n\n");
               exit(1);
          }
          free(irejected);   


          int cap = 10;
          printf ("capped to %d:\t", cap);
          nrejected = fdr(data, cap, alpha, num_tests, &irejected);
          if (nrejected) {
               for (i = 0; i < nrejected; i++) {
                    printf("%f, ", data[irejected[i]]);
               }
          } else {
               printf("None");
          }
          printf ("\n");
                     
          if (nrejected==exp_sig) {/* FIXME should also test values */
               printf("PASS\n\n");
          } else {
               printf("FAIL\n\n");
          }
          free(irejected);   

          printf ("shuffled:\t");
          shuffle(data, NELEMS(data), sizeof(data[0]));
          nrejected = fdr(data, data_len, alpha, num_tests, &irejected);
          if (nrejected) {
               for (i = 0; i < nrejected; i++) {
                    printf("%f, ", data[irejected[i]]);
               }
          } else {
               printf("None");
          }
          printf ("\n");
                     
          if (nrejected==exp_sig) {/* FIXME should also test values */
               printf("PASS\n");
          } else {
               printf("FAIL\n");
          }
          free(irejected);   

          printf ("\n");
     }


     exit(1);

     /* output values according to python implementation 
      *

      # HolmBonferroni
      #
      >>> pvals = [0.010000, 0.010000, 0.030000, 0.050000, 0.005000]
      >>> multiple_testing.HolmBonferroni(pvals).corrected_pvals
      # [0.04, 0.04, 0.06, 0.05, 0.025]      
      >>> multiple_testing.HolmBonferroni(pvals, n=999).corrected_pvals
      # [9.98, 9.98, 29.88, 49.75, 4.995]


      # Bonferroni
      #
      >>> pvals = [0.010000, 0.010000, 0.030000, 0.050000, 0.005000]
      >>> multiple_testing.Bonferroni(pvals).corrected_pvals
      # [0.05, 0.05, 0.15, 0.25, 0.025]
      >>> multiple_testing.Bonferroni(pvals, n=999).corrected_pvals
      # [9.99, 9.99, 29.97, 49.95, 4.995]


      # FDR
      #
      >>> pvals = [0.600000, 0.070000, 0.490000, 0.200000, 0.480000, 0.740000, 0.680000, 0.010000, 0.970000, 0.380000, 0.032000, 0.070000]
      >>> sorted([pvals[i] for i in fdr.fdr(pvals, a=0.20)])
      # [0.01, 0.032]

      */
     for (ntests=-1; ntests<1000; ntests+=1000) {
#if 1
          double data[] = {0.01, 0.01, 0.03, 0.05, 0.005};
          int arrlen = 5;
#else
          double data[] = {0.308799, 0.297089, 0.150936, 0.518048, 0.929977, 0.533344, 0.215166, 0.046900, 0.559717, 0.963137, 0.271862, 0.081730, 0.127197, 0.472981, 0.872270, 0.830436, 0.008807, 0.512721, 0.798971, 0.102211, 0.223481, 0.877118, 0.895051, 0.258633, 0.262466, 0.971096, 0.933613, 0.684905, 0.370404, 0.620535, 0.329913, 0.589691, 0.327275, 0.493765, 0.298086, 0.596776, 0.441615, 0.843791, 0.790259, 0.605214, 0.681068, 0.253347, 0.079298, 0.010991, 0.120785, 0.885083, 0.587810, 0.085836, 0.062631, 0.188596, 0.799452, 0.086936, 0.802441, 0.589041, 0.247903, 0.664755, 0.535578, 0.779242, 0.632830, 0.848579, 0.421346, 0.446816, 0.889726, 0.545608, 0.252037, 0.432924, 0.038197, 0.443023, 0.022724, 0.570956, 0.524829, 0.047716, 0.097964, 0.671701, 0.364718, 0.113121, 0.894922, 0.148989, 0.741134, 0.895151, 0.473917, 0.807252, 0.682185, 0.577869, 0.809614, 0.381577, 0.117793, 0.776623, 0.934201, 0.665000, 0.248839, 0.169519, 0.514693, 0.698509, 0.162058, 0.296180, 0.519535, 0.222471, 0.158160, 0.113808};
          int arrlen = 100;
#endif
          float alpha = 0.05;
          
          printf("Testing holm-bonferroni with ntests=%d alpha=%f\n", ntests, alpha);

          printf ("in: ");
          for (i = 0; i < arrlen; i++) {
               printf("%f, ", data[i]);
          }
          printf ("\n");

          holm_bonf_corr(data, arrlen, alpha, ntests);
          
          printf ("out: ");
          for (i = 0; i < arrlen; i++) {
               printf("%f, ", data[i]);
          }
          printf ("\n\n");
     }

     for (ntests=-1; ntests<1000; ntests+=1000) {
          double data[] = {0.01, 0.01, 0.03, 0.05, 0.005};
          int arrlen = 5;
          /* float alpha = 0.05; */

          printf("Testing bonferroni with ntests=%d...\n", ntests);

          printf ("in: ");
          for (i = 0; i < arrlen; i++) {
               printf("%f, ", data[i]);
          }
          printf ("\n");
          
          bonf_corr(data, arrlen, ntests);
          
          printf ("out: ");
          for (i = 0; i < arrlen; i++) {
               printf("%f, ", data[i]);
          }
          printf ("\n\n");
     }


     for (ntests=-1; ntests<1000; ntests+=1000) {
          /* Test data from : http://udel.edu/~mcdonald/statmultcomp.html
           */          
#if 1
          double data[] = {0.6, 0.07, 0.49, 0.2, 0.48, 0.74,
                            0.68, 0.01, 0.97, 0.38, 0.032, 0.07};
          int arrlen = 12;
          float alpha = 0.20;
#else
          double data[] = {0.308799, 0.297089, 0.150936, 0.518048, 0.929977, 0.533344, 0.215166, 0.046900, 0.559717, 0.963137, 0.271862, 0.081730, 0.127197, 0.472981, 0.872270, 0.830436, 0.008807, 0.512721, 0.798971, 0.102211, 0.223481, 0.877118, 0.895051, 0.258633, 0.262466, 0.971096, 0.933613, 0.684905, 0.370404, 0.620535, 0.329913, 0.589691, 0.327275, 0.493765, 0.298086, 0.596776, 0.441615, 0.843791, 0.790259, 0.605214, 0.681068, 0.253347, 0.079298, 0.010991, 0.120785, 0.885083, 0.587810, 0.085836, 0.062631, 0.188596, 0.799452, 0.086936, 0.802441, 0.589041, 0.247903, 0.664755, 0.535578, 0.779242, 0.632830, 0.848579, 0.421346, 0.446816, 0.889726, 0.545608, 0.252037, 0.432924, 0.038197, 0.443023, 0.022724, 0.570956, 0.524829, 0.047716, 0.097964, 0.671701, 0.364718, 0.113121, 0.894922, 0.148989, 0.741134, 0.895151, 0.473917, 0.807252, 0.682185, 0.577869, 0.809614, 0.381577, 0.117793, 0.776623, 0.934201, 0.665000, 0.248839, 0.169519, 0.514693, 0.698509, 0.162058, 0.296180, 0.519535, 0.222471, 0.158160, 0.113808};
          int arrlen = 100;
          float alpha = 10;
#endif

          int nrejected;
          long int* irejected;/* indices of rejected i.e. significant values */

          printf("Testing fdr with ntests=%d alpha=%f...\n", ntests, alpha);

          printf ("in: ");
          for (i = 0; i < arrlen; i++) {
               printf("%f, ", data[i]);
          }
          printf ("\n");

          nrejected = fdr(data, arrlen, alpha, ntests, &irejected);
          
          printf ("out rejected: ");
          if (nrejected) {
               for (i = 0; i < nrejected; i++) {
                    printf("%f, ", data[irejected[i]]);
               }
          } else {
               printf("None");
          }
          printf ("\n\n");
          
          free(irejected);
     }

     return EXIT_SUCCESS;
}
#endif
