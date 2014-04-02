/* -*- c-file-style: "k&r"; indent-tabs-mode: nil; -*-
 *
 * FIXME Copyright update
 *
 *********************************************************************/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <ctype.h>

/* lofreq includes */
#include "utils.h"
#include "multtest.h"



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


/* gcc -o multtest multtest.c utils.c log.c -ansi -Wall -DMULTTEST_MAIN */
#ifdef MULTTEST_MAIN
int main()
{
     int i;
     int ntests;

     fprintf(stderr, "WARN %s\n", "test holm bonf one last time then upload");

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
          int* irejected;/* indices of rejected i.e. significant values */

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
