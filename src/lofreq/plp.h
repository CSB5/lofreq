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

#ifndef PLP_H
#define PLP_H

#include "htslib/faidx.h"
#include "utils.h"
#include "vcf.h"
#include "utils.h"

/* mpileup configuration flags 
 */
#define MPLP_NO_ORPHAN   0x10
#define MPLP_BAQ         0x20
#define MPLP_REDO_BAQ    0x40
#define MPLP_EXT_BAQ     0x80
#define MPLP_IDAQ        0x100
#define MPLP_REDO_IDAQ   0x200
#define MPLP_USE_SQ      0x400
#define MPLP_ILLUMINA13  0x800


extern const char *bam_nt4_rev_table; /* similar to bam_nt16_rev_table */
#define NUM_NT4 5 /* strlen(bam_nt4_rev_table); */

extern const unsigned char bam_nt4_table[256];


/* mpileup configuration structure 
 */
typedef struct {
     int max_mq, min_mq;
     int flag; /* tag: shared */
     int max_depth;
     int min_plp_bq; /* use with caution: this makes lofreq blind to any bases below this value */
     int min_plp_idq;
     int def_nm_q;
     char *reg;
     char *fa;
     faidx_t *fai;
     void *bed;
     char *alnerrprof_file; /* logically belongs to varcall_conf, but we need it here since only here the bam header is known */
     char cmdline[1024];
} mplp_conf_t;


typedef struct {
     char *target; /* chromsome or sequence name */
     int pos; /* position */
     char ref_base; /* uppercase reference base (given by fasta) */
     char cons_base[MAX_INDELSIZE]; /* uppercase consensus base according to base-counts, after read-level filtering. */
     int coverage_plp; /* original samtools value. upper count limit for all kept values */
     int num_bases; /* number of bases after base filtering */
     /* num_ins and num_dels gives 'num_indels' */
     int num_ign_indels; /* a hack: indels often get filtered because of low quality of missing qualities in bam file. we need to know nevertheless they are present. this is the count of all "ignored" indels */

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

     /* Indel qualities are stored separately according to the type of
      * indel event observed. Insertions and deletions are considered 
      * independently. If there was no indel event observed,
      * the indel quality, indel mapping quality and indel source quality
      * are stored in *_quals, *_map_quals, *_source_quals. Since no
      * indel was observed, there is no indel alignment quality. If
      * an indel event is observed, the qualities are stored in 
      * the hash table to which *_event_counts points to and keyed to the
      * sequence of the indel event. See utils.h for the data structure for storing
      * indel qualities if an indel event is observed. */

     int num_non_indels;/* non-indel events for which we have indel qualities */

     int num_ins, sum_ins;
     int_varray_t ins_quals; 
     int_varray_t ins_map_quals;
     int_varray_t ins_source_quals;
     ins_event *ins_event_counts;

     int num_dels, sum_dels;
     int_varray_t del_quals; 
     int_varray_t del_map_quals;
     int_varray_t del_source_quals;
     del_event *del_event_counts;
     
     /* fw or rv counts for all non-indel events 
      * fw = 0, rv = 1*/
     long int non_ins_fw_rv[2]; 
     long int non_del_fw_rv[2];

     int has_indel_aqs; /* flag, which is only used to make sure that
                           BAM contained alignment quality for indel calls 
                           (all reads with indels should have those).
                           indels are still predicted if missing, but overcalled. */  
     int hrun; /* homopolymer run at (to the right of) current
                * position. if indels are not left aligned and current
                * position is already a homopolymer this will be taken
                * into account. mainly for filtering low af FP indel
                * at the beginning of poly-AT regions. A del GT>G
                * which is in the sequence context of GTTT will
                * receive an hrun value of 3. same for ins G>GT.
                */
     /* changes here should be reflected in plp_col_init, plp_col_free etc. */
} plp_col_t;


#define PLP_COL_ADD_QUAL(p, q)   int_varray_add_value((p), (q))

/* initialize members of preallocated varcall_conf */
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
