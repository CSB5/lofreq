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

#ifndef SNPCALLER_H
#define SNPCALLER_H

#include "vcf.h"
#include "plp.h"
#include "defaults.h"



typedef struct {
     int min_bq;
     int min_alt_bq;
     int def_alt_bq;

     int min_jq;
     int min_alt_jq;
     int def_alt_jq;

     int bonf_dynamic; /* boolean: incr bonf as we go along. eventual
                        * filtering of all has to be done by
                        * caller! */
     int min_cov;
     long long int bonf_subst; /* warning: changed dynamically ! */
     long long int bonf_indel;
     float sig;
     vcf_file_t vcf_out;
     int flag; /* FIXME doc? */

     /* FIXME the following two logically don't belong her but
      * would require a new structure */
     int only_indels; 
     int no_indels; 

} varcall_conf_t;


double
merge_srcq_baseq_and_mapq(const int sq, const int bq, const int mq);

double
merge_srcq_baseq_mapq_and_alnq(const int sq, const int bq, const int mq, const int aq);

void
plp_to_errprobs(double **err_probs, int *num_err_probs, 
                int *alt_bases, int *alt_counts, int *alt_raw_counts,
                const plp_col_t *p, varcall_conf_t *conf);
void 
plp_to_ins_errprobs(double **err_probs, int *num_err_probs, 
                    const plp_col_t *p, varcall_conf_t *conf,
                    char key[MAX_INDELSIZE]);

void 
plp_to_del_errprobs(double **err_probs, int *num_err_probs, 
                    const plp_col_t *p, varcall_conf_t *conf,
                    char key[MAX_INDELSIZE]);

void
init_varcall_conf(varcall_conf_t *c);

void
dump_varcall_conf(const varcall_conf_t *c, FILE *stream) ;


extern double *
poissbin(long double *pvalue, const double *err_probs,
         const int num_err_probs, const int num_failures, 
         const long long int bonf, const double sig);
extern int
snpcaller(long double *snp_pvalues, const double *err_probs,
          const int num_err_probs, const int *noncons_counts,
          const long long int bonf_factor,
          const double sig_level);


#endif
