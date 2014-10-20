/*********************************************************************
 *
 * Copyright (C) 2011-2014 Genome Institute of Singapore
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 *********************************************************************/
#ifndef UTILS_H
#define UTILS_H

#include <limits.h>
#include <float.h>
#include <math.h>
#include <uthash.h>

#define MAX_INDELSIZE 128

#define HAS_GZIP_EXT(f)  (strlen(f)>3 && 0==strncmp(& f[strlen(f)-3], ".gz", 3))


#define PHREDQUAL_TO_PROB(phred) (phred==INT_MAX ? DBL_MIN : pow(10.0, -1.0*(phred)/10.0))

/* requires that prob comes out of our functions is is never zero! */
#define PROB_TO_PHREDQUAL(prob) (int)(-10.0 * log10l(prob))
#define PROB_TO_PHREDQUAL_SAFE(prob) (prob<=0.0 ? INT_MAX : (int)(-10.0 * log10l(prob)))

#define BASECALLQUAL_VALID_RANGE(phred) ((phred)>=0 && (phred)<100)

#define BASENAME(x) strrchr((x), '/') ? strrchr((x), '/')+1 : (x)

int file_exists(const char *fname);
int  is_dir(const char *path);
int ae_load_file_to_memory(const char *filename, char **result);
int int_cmp(const void *a, const void *b);
int dbl_cmp(const void *a, const void *b);
int argmax_d(const double *arr, const int n);
long int count_lines(const char *filename);

typedef struct {
     unsigned long int n; /* number of elements stored */
     int *data; /* actual array of data */

     size_t grow_by_size; /* if needed grow array by this value. will double previous size if <=1 */
     size_t alloced; /* actually allocated size for data */
} int_varray_t;

void int_varray_add_value(int_varray_t *a, const int value);
void int_varray_free(int_varray_t *a);
void int_varray_init(int_varray_t *a, 
                     const size_t grow_by_size);

int
ls_dir(char ***matches, const char *path, const char *pattern,
       const int sort_lexi);

char *
join_paths(char **p1, const char *p2);

void
chomp(char *s);

char *
readlink_malloc(const char *filename);

char *
resolved_path(const char *path);

double
dbl_median(double data[], int size);

int
int_median(int data[], int size);
void
strstrip(char *str);
int
is_newer(const char *p1, const char *p2);

/* utility hash functions for indel calling */

typedef struct {
  char key[MAX_INDELSIZE];
  int count;
  int cons_quals;
  int_varray_t ins_quals;
  int_varray_t ins_aln_quals;
  int_varray_t ins_map_quals;
  int_varray_t ins_source_quals;
  long int fw_rv[2];
  UT_hash_handle hh_ins;
} ins_event;

void add_ins_sequence(ins_event **head_ins_count, char seq[], 
  int ins_qual, int ins_aln_qual, int ins_map_qual, int ins_source_qual, 
  int fw_rv);
ins_event *find_ins_sequence(ins_event *const *head_ins_counts, char seq[]);
void destruct_ins_event_counts(ins_event **head_ins_counts);

typedef struct {
  char key[MAX_INDELSIZE];
  int count;
  int cons_quals;
  int_varray_t del_quals;
  int_varray_t del_aln_quals;
  int_varray_t del_map_quals;
  int_varray_t del_source_quals;
  long int fw_rv[2];
  UT_hash_handle hh_del;
} del_event;

void add_del_sequence(del_event **head_del_counts, char seq[], 
  int del_qual, int del_aln_qual, int del_map_qual, int del_source_qual, 
  int fw_rv);
del_event * find_del_sequence(del_event *const *head_del_counts, char seq[]);
void destruct_del_event_counts(del_event **head_del_counts);

#endif
