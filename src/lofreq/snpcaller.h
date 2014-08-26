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
     int dont_skip_n;
     long long int bonf; /* warning: changed dynamically ! */
     float sig;
     vcf_file_t vcf_out;
     int flag;
} snvcall_conf_t;


double
merge_srcq_baseq_and_mapq(const int sq, const int bq, const int mq);

double
merge_srcq_baseq_mapq_and_alnq(const int sq, const int bq, const int mq, const int aq);

void
plp_to_errprobs(double **err_probs, int *num_err_probs, 
                int *alt_bases, int *alt_counts, int *alt_raw_counts,
                const plp_col_t *p, snvcall_conf_t *conf);

void
init_snvcall_conf(snvcall_conf_t *c);

void
dump_snvcall_conf(const snvcall_conf_t *c, FILE *stream) ;


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
