#ifndef SNPCALLER_H
#define SNPCALLER_H

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

#include "vcf.h"
#include "plp.h"
#include "defaults.h"



typedef struct {
     int min_altbq;
     int def_altbq;
     int bonf_dynamic; /* boolean: incr bonf as we go along. eventual
                        * filtering of all has to be done by
                        * caller! */
     int min_cov;
     int dont_skip_n;
     long long int bonf; /* warning: changed dynamically ! */
     float sig;
     vcf_file_t vcf_out;
     int flag;
} snvcall_conf_t;


void
plp_to_errprobs(double **err_probs, int *num_err_probs, 
                int *alt_bases, int *alt_counts, int *alt_raw_counts,
                const plp_col_t *p, snvcall_conf_t *conf);

void
init_snvcall_conf(snvcall_conf_t *c);

void
dump_snvcall_conf(const snvcall_conf_t *c, FILE *stream) ;


extern double *
poissbin(double *pvalue, const double *err_probs,
         const int num_err_probs, const int num_failures, 
         const long long int bonf, const double sig);
extern int
snpcaller(double *snp_pvalues, const double *err_probs,
          const int num_err_probs, const int *noncons_counts,
          const long long int bonf_factor,
          const double sig_level);


#endif
