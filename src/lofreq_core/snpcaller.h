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

int snpcaller(double *snp_pvalues, const int *phred_quals,
                   const int num_phred_quals, const int *noncons_counts,
                   const unsigned long int bonf_factor,
                   const double sig_level);
