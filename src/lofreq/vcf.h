/* -*- c-file-style: "k&r"; indent-tabs-mode: nil; -*-
 *
 *
 * FIXME missing license
 *
 */
#ifndef VCF_H
#define VCF_H


#include "zlib.h"
#include "uthash.h"


typedef struct {
     FILE *fh;
     gzFile *fh_gz;
     int gz;
} vcf_file_t;

typedef struct {
     char *chrom;
     long int pos; /* zero offset */
     char *id;
     char *ref;
     char *alt;
     int qual;
     char *filter;
     char *info;

     /* genotyping info (not used in lofreq) */
     char *format;
     int num_samples;
     char **samples;
} var_t;

typedef struct {
     int ref_fw;
     int ref_rv;
     int alt_fw;
     int alt_rv;
} dp4_counts_t;


typedef struct {
     char *key; /* according to uthash doc this should be const but then we can't free it */
     var_t *var;
     UT_hash_handle hh;
} var_hash_t;


void var_hash_add(var_hash_t **var_hash, char *key, var_t *var);
void var_hash_free_elem(var_hash_t *hash_elem_ptr);
void  var_hash_free_table(var_hash_t *var_hash);


#define VCF_MISSING_VAL_STR "."
#define VCF_MISSING_VAL_CHAR VCF_MISSING_VAL_STR[0]


#define VCF_VAR_PASSES(v) ((v)->filter[0]==VCF_MISSING_VAL_CHAR || 0==strncmp((v)->filter, "PASS", 4))


#define VCF_EOF(vf) ((vf)->gz ? gzeof((vf)->fh_gz) : feof((vf)->fh))

/* there is no gzvprintf therefore we can't make this a function. this
 * way howver you won't know how many chars were written if any.
 * could add extra variable...
 */
#define VCF_PRINTF(vf, fmt, args...) if (vf->gz) { \
                                       gzprintf(vf->fh_gz, fmt, ## args); \
                                     } else { \
                                       fprintf(vf->fh, fmt, ## args); \
                                     };
int
vcf_file_seek(vcf_file_t *f, long int offset, int whence);
int
vcf_file_open(vcf_file_t *f, const char *path, const int gzip, const char mode);
int
vcf_file_close(vcf_file_t *f);
char *
vcf_file_gets(vcf_file_t *f, int len, char *line);


int vcf_get_dp4(dp4_counts_t *dp4, var_t *var);

void vcf_new_var(var_t **var);
void vcf_free_var(var_t **var);
void vcf_cp_var(var_t **dest, var_t *src);

void vcf_var_key(char **key, var_t *var);
void vcf_var_key_pos_only(char **key, var_t *var);

int vcf_parse_header(char **header, vcf_file_t *vcf_file);
int vcf_skip_header(vcf_file_t *vcf_file);
int vcf_parse_var(vcf_file_t *vcf_file, var_t *var);
int vcf_parse_vars(var_t ***vars, vcf_file_t *vcf_file, int only_passed);

int vcf_var_is_indel(const var_t *var);
int vcf_var_has_info_key(char **value, const var_t *var, const char *key);
int vcf_var_filtered(const var_t *var);
char *vcf_var_add_to_filter(var_t *var, const char *filter_name);
char *vcf_var_add_to_info(var_t *var, const char *info_str);
void vcf_var_sprintf_info(var_t *var,
                          const int *dp, const float *af, const int *sb,
                          const dp4_counts_t *dp4,
                          const int is_indel, const int is_consvar);
void vcf_write_var(vcf_file_t *vcf_file, const var_t *var);
void vcf_write_header(vcf_file_t *vcf_file, const char *header);
void vcf_write_new_header(vcf_file_t *vcf_file, const char *srcprog, const char *reffa);
void vcf_header_add(char **header, const char *info);
#endif
