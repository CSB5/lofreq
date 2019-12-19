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

#ifndef SAMUTILS_H
#define SAMUTILS_H

#include "htslib/sam.h"




typedef enum {
        OP_MATCH,
        OP_MISMATCH,
        OP_INS,
        OP_DEL,
        NUM_OP_CATS,
} op_cat_t;

#define STR(name) # name

static char *op_cat_str[] = {
    STR(OP_MATCH),
    STR(OP_MISMATCH),
    STR(OP_INS),
    STR(OP_DEL),
    STR(NUM_OP_CATS)
};


char *
cigar_str_from_bam(const bam1_t *b);

int
count_cigar_ops(int *counts, int **quals,
                const bam1_t *b, const char *ref, int min_bq,
                char *target);


#ifdef USE_ALNERRPROF

typedef struct {
     int num_targets; /* bam_header->n_targets */
     int *prop_len; /* one prop length per target: index is tid */
     double **props; /* one prop array per target: index is tid */
} alnerrprof_t;


void
normalize_alnerrprof(alnerrprof_t *alnerrprof);

int
parse_alnerrprof_statsfile(alnerrprof_t *alnerrprof, const char *path, bam_hdr_t *bam_header);

void
calc_read_alnerrprof(double *alnerrprof, unsigned long int *used_pos, 
                        const bam1_t *b, const char *ref);

void
write_alnerrprof_stats(char *target_name, unsigned long int *alnerrprof_usedpos, 
                    double *alnerrprof, int max_obs_read_len, FILE *out);

void
free_alnerrprof(alnerrprof_t *alnerrprof);

#endif

int checkref(char *fasta_file, char *bam_file);

#endif
