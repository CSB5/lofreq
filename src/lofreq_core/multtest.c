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
     int i;
} ixp; /* indexed p-value */

int
ixp_dbl_cmp(const void *a, const void *b)
{
     const ixp *ia = (const ixp *)a;
     const ixp *ib = (const ixp *)b;
     return dbl_cmp(&ia->p, &ib->p);
}


void
bonf_corr(double array[], int asize, int num_tests, double alpha)
{
/* FIXME */
}

void
holm_bonf_corr(double array[], int n, double alpha)
{

     /* FIXME add num_tests appropriately */

     ixp iarr[n];
     int i;
     /* first index the pvalues and store them in the ixp array */
     for (i = 0; i < n; i++) {
          iarr[i].i = i;
          iarr[i].p = array[i];
     }
     qsort(iarr, n, sizeof(ixp), ixp_dbl_cmp);

     int lp = n;
     double pp = iarr[0].p;
     double tp;
     for (i = 0; i < n; i++) {
          /* if the pvalue is different
           * update lp according to how many pvalues
           * have been seen and update the
           * previous pvalue seen */
          if (dbl_cmp(&iarr[i].p, &pp) != 0) {
               lp = n - i;
               pp = iarr[i].p;
          }
          /* if below alpha, correct the
           * pvalue in the original array */
          tp = iarr[i].p * 1./ lp;
          if (dbl_cmp(&tp, &alpha) < 0) {
               array[iarr[i].i] = iarr[i].p * lp;
          }
     }
}

int
fdr(double array[], int n, double alpha, int **irejected)
{

     /* FIXME add num_tests appropriately */

     ixp iarr[n];
     int i;
     /* first index the pvalues and store them in ixp array to sort indices*/
     for (i = 0; i < n; i++) {
          iarr[i].i = i;
          iarr[i].p = array[i];
     }
     qsort(iarr, n, sizeof(ixp), ixp_dbl_cmp);

     int nrejected;
     /* starting from the largest rank, evaluate p(m) < alpha * (m/M) where m is the
      * rank and M is the total number of pvalues. If true, reject pvals 1..m */
     for (i = n; i > 0; i--) { // ranks are 1-based
          if (iarr[i-1].p < (alpha*i/n)) {
               nrejected = i; // therefore, nrejected includes first rejected (0-based)
               break;
          }
     }

     /* return array of indices to rejected pvalues */
     /* FIXME: warning: ‘nrejected’ is used uninitialized in this function */
     (*irejected) = (int*) malloc(nrejected*sizeof(int));
     for (i = 0; i < nrejected; i++) {
          // printf("%d\t%f\t%d\n", iarr[i].i, iarr[i].p);
          (*irejected)[i] = iarr[i].i;
     }
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


#ifdef MULTTEST_MAIN
int main()
{

     printf("Testing holm-bonferroni...\n");

     double array[] = {0.01, 0.01, 0.03, 0.05, 0.005};
     int arrlen = 5;
     float alpha = 0.05;

     holm_bonf_corr(array, arrlen, alpha);

     int i;
     for (i = 0; i < arrlen; i++) {
          printf("%f\n", array[i]);
     }

     /* Test data from : http://udel.edu/~mcdonald/statmultcomp.html
      */

     printf("Testing fdr...\n");
     double array2[] = {0.6, 0.07, 0.49, 0.2, 0.48, 0.74,
                       0.68, 0.01, 0.97, 0.38, 0.032, 0.07};
     int arrlen2 = 12;
     float alpha2 = 0.20;

     int* irejected;
     int nrejected = fdr(array2, arrlen2, alpha2, &irejected);


     for (i = 0; i < nrejected; i++) {
          printf("%f\n", array2[irejected[i]]);
     }

     free(irejected); /* free me! */

}
#endif
