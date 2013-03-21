/* -*- mode: c; tab-width: 4; c-basic-offset: 4;  indent-tabs-mode: nil -*- */

#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <limits.h>

#include "log.h"
#include "utils.h"

#ifndef SIZE_MAX
#define SIZE_MAX ((size_t)-1)
#endif

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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 *********************************************************************/



/* overflow safe int comparison for e.g. qsort.
 *
 * a simple return *ia - *ib; can in theory overflow see
 * http://stackoverflow.com/questions/6103636/c-qsort-not-working-correctly
 */
int int_cmp(const void *a, const void *b)
{
     const int ia = *(const int *)a;
     const int ib = *(const int *)b;
     return ia<ib ? -1 : ia>ib? 1 : 0;
}


/* return index for max double in array. will return the lower index
 * on tie */
int argmax_d(const double *arr, const int n)
{
  int i;
  int maxidx = 0;

  for (i=0; i<n; i++) {
       if (arr[i] > arr[maxidx]) {
            maxidx = i;
       }
  }
  return maxidx;
}


void int_varray_free(int_varray_t *a) 
{
    assert(NULL != a);

    free(a->data); /* save even if a->data==NULL */
    a->data = NULL;
    a->n = a->alloced = a->grow_by_size = 0;
}

void int_varray_init(int_varray_t *a, 
                     const size_t grow_by_size)
{
    assert(NULL != a);

    a->n = 0;
    a->data = NULL;
    a->grow_by_size = grow_by_size;
    a->alloced = 0;
}

void int_varray_add_value(int_varray_t *a, const int value)
{
    assert(NULL != a);

    if (a->n * sizeof(int) == a->alloced) {
        size_t size_to_alloc;        
        if (1 >=  a->grow_by_size) {
             assert(SIZE_MAX - a->alloced > a->alloced);
             size_to_alloc = 0==a->n ? sizeof(int) : a->alloced*2;
        } else {
             assert(SIZE_MAX - a->alloced > a->grow_by_size);
             size_to_alloc = a->alloced + a->grow_by_size;
        }
        a->data = realloc(a->data, size_to_alloc);
        a->alloced = size_to_alloc;
    }
    a->data[a->n] = value;
    a->n++;
}



int file_exists(const char *fname) 
{
  /* from 
   * http://stackoverflow.com/questions/230062/whats-the-best-way-to-check-if-a-file-exists-in-c-cross-platform 
   */
  if (access(fname, F_OK) != -1) {
      return 1;
  } else {
      return 0;
  }
}

/* from http://www.anyexample.com/programming/c/how_to_load_file_into_memory_using_plain_ansi_c_language.xml
 *
 * returns file size (number of bytes) on success or negative number
 * on error
 * 
 * warnings:
 * Function ae_load_file_to_memory returns loaded data size which does not take into account last null-terminate symbol.
 * If you want to use this function to process string data, note that it may work incorrectly with multibyte encodings.
 */
int ae_load_file_to_memory(const char *filename, char **result) 
{ 
	int size = 0;
	FILE *f = fopen(filename, "rb");
	if (f == NULL) { 
		*result = NULL;
		return -1; /* -1 means file opening fail */
	} 
	fseek(f, 0, SEEK_END);
	size = ftell(f);
	fseek(f, 0, SEEK_SET);
	*result = (char *)malloc(size+1);
    if (NULL==result) {
        return -2;
    }
	if (size != fread(*result, sizeof(char), size, f)) {
		free(*result);
		return -3; /* -2 means file reading fail  */
	} 
	fclose(f);
	(*result)[size] = 0;
	return size;
}

/* count number of lines advise from
 * http://stackoverflow.com/questions/8689344/portable-end-of-line-newline-in-c:
 * open in binary mode and count \n. in text mode \n is replaced by a
 * platform specific ELS.
 *
 * Returns value <0 on failure. Otherwise line count.
 */
long int
count_lines(const char *filename)
{
    int c;
    long int count = 0;
	FILE *f = fopen(filename, "rb");

	if (f == NULL) { 
		return -1; /* -1 means file opening fail */
	}
    while (EOF != (c=getc(f))) {
        if ('\n'==c) {
            if (count==LONG_MAX) {
                LOG_FATAL("%s\n", "count overflow!");
                return -2;
            }
            count++;
        }
    }
    fclose(f);
    return count;
}
/* count_lines */


/* returns -1 on error 
 * FIXME we should be using libbam's bed_read()
 */
long long int
bed_pos_sum(const char *bed_file) {
#define BUF_SIZE 1024
    long long int sum = 0;
    char line[BUF_SIZE];
    FILE *fh;

    fh = fopen(bed_file, "r");
    if (NULL == fh) {
        LOG_ERROR("Couldn't open bed-file %s\n", bed_file);
        return -1;
    }

    while (NULL != fgets(line, BUF_SIZE, fh)) {
        char chrom[BUF_SIZE];
        long int start, end;
        if (line[0]=='#') {
            continue;
        }

        if (1 == strlen(line)) {
            LOG_WARN("Skippping empty line in bed-file %s", bed_file);
            continue;
        }

        /* this works with any number of tabs and white-spaces */
        if (3 != sscanf(line, "%s %ld %ld \n", chrom, &start, &end)) {
            LOG_FATAL("Couldn't parse the following line"
                      " from bed-file %s: %s", bed_file, line);
            fclose(fh);
            return -1;
        }
        if (start>end) {
            LOG_FATAL("start > end in the following line"
                      " from bed-file %s: %s", bed_file, line);
            fclose(fh);
            return -1;
        }
        if (sum > LLONG_MAX - (end-start)) {
            LOG_FATAL("%s\n", "count overflow!");
            return -1;
        }
        sum += (end-start);
    }
    fclose(fh);
    return sum;
}
/* bed_pos_sum */

