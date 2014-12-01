/* -*- mode: c; tab-width: 4; c-basic-offset: 4;  indent-tabs-mode: nil -*- */
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



#include <stdlib.h>
#include <stdio.h>

#include "cdflib.h"
#include "binom.h"



#define MIN(X,Y) ((X) < (Y) ? (X) : (Y))



/**
 * @brief Compute cdf and sf
 *
 * P is the cdf evaluated at X, Q is the compliment of the cdf
 * evaluated at X, i.e. 1-P (AKA sf)
 *
 * Returns non-zero status on failure
 *
 */
int binom(double *p, double *q,
          int num_trials, int num_success, double prob_success) 
{
		int which=1;
		int status=1; /* error by default */
		double ompr = 1.0 - prob_success;
		double bound;
        double q2, p2;

		double s = (double)num_success;
		double xn = (double)num_trials;
		double pr = (double)prob_success;

        /* P is always the cdf evaluated at X, Q is always the compliment of the
           cdf evaluated at X, i.e.  1-P, and X is always the value at which the
           cdf  is evaluated. */

		(void) cdfbin(&which, p?p:&p2, q?q:&q2,
			   &s, &xn, &pr, &ompr,
			   &status, &bound);
        
#ifdef DEBUG

		fprintf(stderr, "DEBUG(%s:%s:%d): in num_success = %d\n", 
                __FILE__, __FUNCTION__, __LINE__, num_success);
		fprintf(stderr, "DEBUG(%s:%s:%d): in num_trials = %d\n", 
                __FILE__, __FUNCTION__, __LINE__, num_trials);
		fprintf(stderr, "DEBUG(%s:%s:%d): in pr = %g\n", 
                __FILE__, __FUNCTION__, __LINE__, prob_success);
		fprintf(stderr, "DEBUG(%s:%s:%d): out p=%g\n", 
                __FILE__, __FUNCTION__, __LINE__, p?*p:p2);
		fprintf(stderr, "DEBUG(%s:%s:%d): out q=%g\n", 
                __FILE__, __FUNCTION__, __LINE__, q?*q:q2);
		fprintf(stderr, "DEBUG(%s:%s:%d): out status=%d\n",
                __FILE__, __FUNCTION__, __LINE__, status);
		fprintf(stderr, "DEBUG(%s:%s:%d): out bound=%g\n", 
                __FILE__, __FUNCTION__, __LINE__, bound);
#endif
		
		return status;
}
/* end of binom */



#ifdef BINOM_MAIN


/* 
gcc -pedantic -Wall -g -std=gnu99 -O2 -DBINOM_MAIN -I../cdflib90/ -o binom binom.c utils.c log.c ../cdflib90/libcdf.a -lm
*/
#include <stdlib.h>
#include "log.h"

int main(int argc, char *argv[]) {
     int num_success;
     int num_trials;
     double prob_success;
     double sf_pvalue;
     double cdf_pvalue;

     if (argc<4) {
         fprintf(stderr, "need num_success num_trials and prob_success as args");
         return -1;
     }

     num_success = atoi(argv[1]);
     num_trials = atoi(argv[2]);
     prob_success = atof(argv[3]);


     fprintf(stdout, "num_success=%d num_trials=%d prob_success=%f\n", num_success, num_trials, prob_success);
     if (0 != binom(&cdf_pvalue, &sf_pvalue, num_trials, num_success, prob_success)) {
         fprintf(stderr, "%s\n", "binom() failed");
         return EXIT_FAILURE;
     }
     
     printf("sf: %g\tcdf: %g\n", sf_pvalue, cdf_pvalue);

     printf("sf should be identical to scipy.stats.binom.sf(%d, %d, %f)\n", num_success, num_trials, prob_success);
     printf("cdf should be identical to scipy.stats.binom.cdf(%d, %d, %f)\n", num_success, num_trials, prob_success);
     return EXIT_SUCCESS;
}
#endif
