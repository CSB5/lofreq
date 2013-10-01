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

extern double *
poissbin(double *pvalue, const double *err_probs,
         const int num_err_probs, const int num_failures, 
         const long long int bonf, const double sig);
extern int
snpcaller(double *snp_pvalues, const double *err_probs,
          const int num_err_probs, const int *noncons_counts,
          const long long int bonf_factor,
          const double sig_level);

extern int
source_qual(const bam1_t *b, const char *ref);

#endif
