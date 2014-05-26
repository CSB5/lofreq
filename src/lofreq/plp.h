/* -*- c-file-style: "k&r"; indent-tabs-mode: nil; -*- */
#ifndef PLP_H
#define PLP_H

#include "faidx.h"
#include "utils.h"
#include "vcf.h"


/* mpileup configuration flags 
 */
#define MPLP_NO_ORPHAN   0x10
#define MPLP_REALN       0x20
#define MPLP_REDO_BAQ    0x40
#define MPLP_ILLUMINA13  0x80
#define MPLP_IGNORE_RG   0x100
#define MPLP_USE_SQ      0x200


extern const char *bam_nt4_rev_table; /* similar to bam_nt16_rev_table */
#define NUM_NT4 5 /* strlen(bam_nt4_rev_table); */

extern const unsigned char bam_nt4_table[256];


/* mpileup configuration structure 
 */
typedef struct {
     int max_mq, min_mq;
     int flag; /* tag: shared */
     int capQ_thres;
     int max_depth;
     int min_bq;
     int def_nm_q;
     char *reg;
     char *fa;
     faidx_t *fai;
     void *bed;
     char *alnerrprof_file; /* logically belongs to snvcall_conf, but we need it here since only here the bam header is known */
     char cmdline[1024];
} mplp_conf_t;


typedef struct {
     char *target; /* chromsome or sequence name */
     int pos; /* position */
     char ref_base; /* uppercase reference base (given by fasta) */
     char cons_base; /* uppercase consensus base according to base-counts, after read-level filtering. */
     int coverage; /* coverage after read-level filtering i.e. same as in samtools mpileup (n_plp) but without indels! */

     /* list of qualities: keeping them all here in one place so that
      * filtering can become separate step. alternative is to filter
      * during pileup. the latter doesn't work if you want to filter
      * based on a consensus which you don't know in advance */
     int_varray_t base_quals[NUM_NT4]; 
     int_varray_t baq_quals[NUM_NT4]; 
     int_varray_t map_quals[NUM_NT4]; 
     int_varray_t source_quals[NUM_NT4]; 
#ifdef USE_ALNERRPROF
     int_varray_t alnerr_qual[NUM_NT4]; /* FIXME this should be precomputed and then build into model */
#endif
     long int fw_counts[NUM_NT4]; 
     long int rv_counts[NUM_NT4]; 
     /* fw_counts[b] + rv_counts[b] = x_quals.n = coverage */

     int num_heads; /* number of read starts at this pos */
     int num_tails; /* number of read ends at this pos */

     int num_ins, sum_ins;
     int_varray_t ins_quals; 

     int num_dels, sum_dels;
     int_varray_t del_quals; 

     /* changes here should be reflected in plp_col_init, plp_col_free etc. */
} plp_col_t;


#define PLP_COL_ADD_QUAL(p, q)   int_varray_add_value((p), (q))

/* initialize members of preallocated snvcall_conf */
void init_mplp_conf(mplp_conf_t *c);

int
base_count(const plp_col_t *p, char base);

void
dump_mplp_conf(const mplp_conf_t *c, FILE *stream);

int
mpileup(const mplp_conf_t *mplp_conf, 
        void (*plp_proc_func)(const plp_col_t*, void*),
        void *plp_proc_conf, 
        const int n, const char **fn);

int
source_qual_load_ign_vcf(const char *vcf_path, void *bed);

void
source_qual_free_ign_vars();

int 
var_in_ign_list(var_t *var);

#endif
