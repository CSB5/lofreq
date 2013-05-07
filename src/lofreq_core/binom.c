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

#include "cdflib.h"
#include "binom.h"

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



#ifdef BINOM_MAIN


/* 
gcc -pedantic -Wall -g -std=gnu99 -O2 -DBINOM_MAIN -I../cdflib90/ -o binom binom.c utils.c log.c ../cdflib90/libcdf.a
*/
#include <stdlib.h>
#include "log.h"

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
     if (0 != binom_sf(&pvalue, num_trials, num_success, succ_prob)) {
         LOG_ERROR("%s\n", "binom_sf() failed");
         return EXIT_FAILURE;
     }

     printf("%g\n", pvalue);

     printf("Should be identical to scipy.stats.binom.sf(%d, %d, %f)\n", num_success, num_trials, succ_prob);
     return EXIT_SUCCESS;
}
#endif
