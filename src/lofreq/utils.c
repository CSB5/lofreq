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


#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <limits.h>
#include <dirent.h>
#include <errno.h>
#include <sys/stat.h>
#include <libgen.h>
#include <ctype.h>

#include "log.h"
#include "utils.h"

#ifndef SIZE_MAX
#define SIZE_MAX ((size_t)-1)
#endif

#define DIR_SEP "/"



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


int dbl_cmp(const void *a, const void *b)
{
     const double da = *(const double *)a;
     const double db = *(const double *)b;

     /* epsilon stuff needed/working at all? */
     if (fabs(da-db) < DBL_EPSILON) {
          return 0;
     }
     return da<db ? -1 : da>db? 1 : 0;
}


int str_cmp(const void *a, const void *b)
{ 
    return strcmp(*(char **)a, *(char **)b);
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


/* returns 1 of path is directory, otherwise 0 if it's anything else
 * or if there's permission problem 
*/
int is_dir(const char *path)
{
     struct stat s;
     if (stat(path, &s) == 0) {
          if (s.st_mode & S_IFDIR) {
               return 1;
          } else {
               return 0;
          }
     } else {
          /* error: could check set errno for more details */
          return 0;
     }         

}

/* also exists in htslib. see http://en.wikipedia.org/wiki/Weak_symbol */
#pragma weak file_exists
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


/* returns -1 on error, otherwise number of matches. caller has to
 * free matches */
int
ls_dir(char ***matches, const char *path, const char *pattern,
       const int sort_lexi)
{
    DIR* d = opendir(path);
    struct dirent *sd = NULL;
    int num_matches = 0;

    (*matches) = NULL;

    if (d == NULL) {
        LOG_ERROR("Couldn't open path %s\n", path);
        return -1;
    }

    while (NULL != (sd = readdir(d))) {/* readdir not thread safe */
        int match = 0;
        if (pattern && strstr(sd->d_name, pattern)) {
            match = 1;
        } else if (NULL==pattern) {
            match = 1;
        }
        if (0 == match) {
            continue;
        }
        num_matches += 1;

        (*matches) = realloc((*matches), num_matches*sizeof(char*)); /* FIXME inefficient one by one allocation */
        if (NULL == (*matches)) {
            LOG_ERROR("%s\n", "Realloc failed");
            return -1;
        }
        (*matches)[num_matches-1] = calloc(strlen(path) +
                                        strlen(sd->d_name) +
                                        1 /*/*/ +1 /*\0*/,
                                        sizeof(char));
        sprintf((*matches)[num_matches-1], "%s/%s", path, sd->d_name);
    }
    closedir(d);

    if (sort_lexi) {
        qsort((*matches), num_matches, sizeof(char*), *str_cmp);
    }
    return num_matches;
}



/* appends dir p2 to p1 and canonicalizes the pathname. returns NULL
 * on error or if normalized path doesn't exist. will allocate memory
 * for p1 as needed.
 */
char * 
join_paths(char **p1, const char *p2) {
     int bufsize;
     char *buf;
     char *buf_resolved;

     if (NULL == p1 || NULL == p2) {
          return NULL;
     }

     bufsize = strlen(*p1) + 1 + strlen(p2) + 1;
     if (bufsize < PATH_MAX) {
          bufsize = PATH_MAX; /* realpath requirement */
     }
     buf = malloc(bufsize * sizeof(char));
     buf_resolved = malloc(bufsize * sizeof(char));

     buf[0] = '\0';
     (void) strcat(buf, *p1);
     (void) strcat(buf, DIR_SEP);
     (void) strcat(buf, p2);
     if (NULL == realpath(buf, buf_resolved)) {
#if 0
          LOG_WARN("Couldn't normalize %s: %s\n",
                   buf, strerror(errno));
#endif
          free(buf_resolved);
          free(buf);
          return NULL;
     } 
     *p1 = realloc(*p1, (strlen(buf_resolved)+1)*sizeof(char));
     (void) strcpy(*p1, buf_resolved);

     free(buf_resolved);
     free(buf);

     return *p1;
}



/* taken from
 * http://www.delorie.com/gnu/docs/glibc/libc_279.html
 * needed because if readlink's 'return value equals size, you cannot
 * tell whether or not there was room to return the entire name'.
 */
char *
readlink_malloc(const char *filename)
{
     int size = 100;
     char *buffer = NULL;
     
     while (1) {
          int nchars = readlink(filename, buffer, size);
          buffer = (char *)realloc(buffer, size);
          if (nchars < 0) {
               free(buffer);
               return NULL;
          }
          if (nchars < size) {
               return buffer;
          }
          size *= 2;
     }
}


/* follows symlinks until resolved and returns realpath. returns NULL
 * on error, otherwise true path. caller has to free
 */
char *
resolved_path(const char *path)
{
     char *resolved_path, *tmp_path;
     char orig_wd[PATH_MAX];

     if (NULL == getcwd(orig_wd, PATH_MAX)) {
          return NULL;
     }

     resolved_path = strdup(path);
     while (1) {
          char realpath_buf[PATH_MAX];
          struct stat stat_buf;

          if (lstat(resolved_path, &stat_buf)) {
               /*LOG_WARN("%s\n", "lstat() failed");*/
               free(resolved_path);
               resolved_path = NULL;
               goto chdir_and_return;
          }

          /* done if not a link */
          if (! S_ISLNK(stat_buf.st_mode)) {
               /*LOG_FIXME("no more link: %s\n", resolved_path);*/
               break;
          }

          /* read link and change to dirname of link */
          if (NULL == (tmp_path = readlink_malloc(resolved_path))) {
               LOG_ERROR("%s\n", "readlink() failed.");
               free(resolved_path);
               resolved_path = NULL;
               goto chdir_and_return;
          }
          if (-1 == chdir(dirname(resolved_path))) {
               LOG_ERROR("%s\n", "chdir() failed.");
               free(tmp_path);
               free(resolved_path);
               return NULL;
          }
          /*LOG_FIXME("Now in %s\n", dirname(resolved_path));*/

          if (NULL == realpath(tmp_path, realpath_buf)) {
               LOG_ERROR("realpath failed on %s\n", tmp_path);
               free(tmp_path);
               free(resolved_path);
               resolved_path = NULL;
               goto chdir_and_return;
          }
                    
          free(tmp_path);
          free(resolved_path);
          resolved_path = strdup(realpath_buf);
     }

chdir_and_return:

     if (-1 == chdir(orig_wd)) {
          LOG_ERROR("%s\n", "chdir() failed. Trying to continue...");
     }
     /*LOG_FIXME("resolved_path is now %s\n", resolved_path);*/

     return resolved_path;
}

/* FIXME use wirth's method instead for larger arrays 
 * FIXME Make malloc optional in case input data can be sorted
 */
int
int_median(int data[], int size)
{
     int ret;
     int *sdata;

     if (size==0) {
          return 0;
     }
     sdata = malloc(sizeof(int) * size);
     memcpy(sdata, data, sizeof(int) * size);
     qsort(sdata, size, sizeof(int), int_cmp);
     if (size%2 == 0) {
          /* even number: return mean of the two elements in the middle */
          ret = (sdata[size/2] + sdata[size/2 - 1]) / 2.0;
     } else {
          /* odd number: return element in middle */
          ret = sdata[size/2];
     }

     free(sdata);
     return ret;
}



/* FIXME use wirth's method instead for larger arrays
 * FIXME Make malloc optional in case input data can be sorted */
double
dbl_median(double data[], int size)
{
     double ret;
     double *sdata;

     if (size==0) {
          return 0.0;
     }
     sdata =  malloc(sizeof(double) * size);
     memcpy(sdata, data, sizeof(double) * size);
     qsort(sdata, size, sizeof(double), dbl_cmp);
     if (size%2 == 0) {
          /* even number: return mean of the two elements in the middle */
          ret = (sdata[size/2] + sdata[size/2 - 1]) / 2.0;
     } else {
          /* odd number: return element in middle */
          ret = sdata[size/2];
     }
     free(sdata);
     return ret;
}


void chomp(char *s)
{
     if (!s) {
          return;
     }
     int end = strlen(s)-1;
     while (end >= 0 && (s[end]=='\n' || s[end]=='\r')) {
          s[end]='\0';
          end = end-1;
     }
}



void
strstrip(char *str)
{
     size_t size;

     fprintf(stderr, "FIXME untested function\n"); exit(1);
     if (! str) {
          return;
     }

     size = strlen(str);
     if (!size) {
          return;
     }

     /* rstrip */
     while (size>0 && isspace(str[size-1])) {
          str[--size] = 0;
     }
     /* lstrip */
     while (*str && isspace(*str)) {
          str++;
     }
}


/* check if first file is newer (mtime) than second. returns 1 if yes,
 * 0 if not and -1 on error */
int
is_newer(const char *p1, const char *p2)
{
     struct stat s1;
     struct stat s2;
     int res;

     if (!p1 || !p2) {
          return -1;
     }
     res = stat(p1, &s1);
     if (res == 0 && !S_ISREG(s1.st_mode)) {
          /* exists but not regular file */
          return -1;
     } else if (res != 0) {
          return -1;
     }

     res = stat(p2, &s2);
     if (res == 0 && !S_ISREG(s2.st_mode)) {
          /* exists but not regular file */
          return -1;
     } else if (res != 0) {
          /* stat failed */
          return -1;
     }

     return s1.st_mtime > s2.st_mtime;
}

void add_ins_sequence(ins_event **head_ins_counts, char seq[], 
     int ins_qual, int ins_aln_qual, int ins_map_qual, int ins_source_qual, 
     int fw_rv) {
     ins_event *it = NULL;
     int seq_length = strlen(seq);
     const int grow_by_size = 16384;
    
     HASH_FIND(hh_ins, *head_ins_counts, seq, seq_length, it);
     if (it) {
          it->count += 1;
          it->cons_quals += ins_qual;
          
          it->fw_rv[fw_rv] += 1;
          
          int_varray_add_value(& it->ins_quals, ins_qual);
          int_varray_add_value(& it->ins_aln_quals, ins_aln_qual);
          int_varray_add_value(& it->ins_map_quals, ins_map_qual);
          int_varray_add_value(& it->ins_source_quals, ins_source_qual);

     } else {
          it = malloc(sizeof(ins_event));
          strncpy((char *)it->key, seq, MAX_INDELSIZE-1);
          it->count = 1;
          it->cons_quals = ins_qual;
          
          it->fw_rv[0] = it->fw_rv[1] = 0;
          it->fw_rv[fw_rv] += 1;

          int_varray_init(& it->ins_quals, grow_by_size);
          int_varray_init(& it->ins_aln_quals, grow_by_size);
          int_varray_init(& it->ins_map_quals, grow_by_size);
          int_varray_init(& it->ins_source_quals, grow_by_size);
          
          int_varray_add_value(& it->ins_quals, ins_qual);
          int_varray_add_value(& it->ins_aln_quals, ins_aln_qual);
          int_varray_add_value(& it->ins_map_quals, ins_map_qual);
          int_varray_add_value(& it->ins_source_quals, ins_source_qual);

          HASH_ADD_KEYPTR(hh_ins, *head_ins_counts, it->key, seq_length, it);
     }
}

ins_event * find_ins_sequence(ins_event *const *head_ins_counts, char seq[]) {
     ins_event *it = NULL;
     HASH_FIND(hh_ins, *head_ins_counts, seq, strlen(seq), it);
     return it;
}

void destruct_ins_event_counts(ins_event **head_ins_counts) {
     ins_event *it_ins, *it_tmp;
     HASH_ITER(hh_ins, *head_ins_counts, it_ins, it_tmp) {
          HASH_DELETE(hh_ins, *head_ins_counts, it_ins);
          int_varray_free(& it_ins->ins_quals);
          int_varray_free(& it_ins->ins_aln_quals);
          int_varray_free(& it_ins->ins_map_quals);
          int_varray_free(& it_ins->ins_source_quals);
          free(it_ins);
     }
}

void add_del_sequence(del_event **head_del_counts, char seq[], 
     int del_qual, int del_aln_qual, int del_map_qual, int del_source_qual, 
     int fw_rv) {
     del_event *it = NULL;
     int seq_length = strlen(seq);
     const int grow_by_size = 16384;

     HASH_FIND(hh_del, *head_del_counts, seq, seq_length, it);
     if (it) {
          it->count += 1;
          it->cons_quals += del_qual;
          
          it->fw_rv[fw_rv] += 1;
          
          int_varray_add_value(& it->del_quals, del_qual);
          int_varray_add_value(& it->del_aln_quals, del_aln_qual);
          int_varray_add_value(& it->del_map_quals, del_map_qual);
          int_varray_add_value(& it->del_source_quals, del_source_qual);
     
     } else {
          it = malloc(sizeof(del_event));
          strncpy((char *)it->key, seq, MAX_INDELSIZE-1);
          it->count = 1;
          it->cons_quals = del_qual;
          
          it->fw_rv[0] = it->fw_rv[1] = 0;
          it->fw_rv[fw_rv] += 1;

          int_varray_init(& it->del_quals, grow_by_size);
          int_varray_init(& it->del_aln_quals, grow_by_size);
          int_varray_init(& it->del_map_quals, grow_by_size);
          int_varray_init(& it->del_source_quals, grow_by_size);

          int_varray_add_value(& it->del_quals, del_qual);
          int_varray_add_value(& it->del_aln_quals, del_aln_qual);
          int_varray_add_value(& it->del_map_quals, del_map_qual);
          int_varray_add_value(& it->del_source_quals, del_source_qual);

          HASH_ADD_KEYPTR(hh_del, *head_del_counts, it->key, seq_length, it);
     }
}

del_event * find_del_sequence(del_event *const *head_del_counts, char seq[]) {
     del_event *it = NULL;
     HASH_FIND(hh_del, *head_del_counts, seq, strlen(seq), it);
     return it;
}

void destruct_del_event_counts(del_event **head_del_counts) {
     del_event *it_del, *it_tmp;
     HASH_ITER(hh_del, *head_del_counts, it_del, it_tmp) {
          HASH_DELETE(hh_del, *head_del_counts, it_del);
          int_varray_free(& it_del->del_quals);
          int_varray_free(& it_del->del_aln_quals);
          int_varray_free(& it_del->del_map_quals);
          int_varray_free(& it_del->del_source_quals);
          free(it_del);
     }
}

void strtoupper(char *s) {
     for (; *s != '\0'; s++) {
          *s = toupper(*s);
     }
}


/* gcc -o utils utils.c log.c -DXMAIN -Wall -ansi -pedantic  */
#ifdef MEDIAN_MAIN
int main(int argc, char **argv)
{
     int i;
     double *data;

     data = malloc((argc-1) * sizeof(double));
     for (i=1; i<argc; i++) {
          printf("%f\n", atof(argv[i]));
          data[i-1] = atof(argv[i]);
     }
     printf("median = %f\n", dbl_median(data, argc-1));
     free(data);
     return 0;
}
#endif

#ifdef  NEWER_MAIN
int main(int argc, char **argv)
{
    printf("is_newer %s %s = %d\n", argv[1], argv[2], is_newer(argv[1], argv[2]));
     return 0;
}
#endif
