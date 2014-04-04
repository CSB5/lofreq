#ifndef UTILS_H
#define UTILS_H

#include <limits.h>
#include <float.h>
#include <math.h>

#define HAS_GZIP_EXT(f)  (strlen(f)>3 && 0==strncmp(& f[strlen(f)-3], ".gz", 3))


#define PHREDQUAL_TO_PROB(phred) (phred==INT_MAX ? DBL_MIN : pow(10.0, -1.0*(phred)/10.0))

#if 0
#define PROB_TO_PHREDQUAL(prob) (prob<0.0+DBL_EPSILON ? INT_MAX : (int)(-10.0 * log10(prob)))
#endif
/* requires that prob comes out of our functions is is never zero! */
#define PROB_TO_PHREDQUAL(prob) (int)(-10.0 * log10(prob))
#define PROB_TO_PHREDQUAL_SAFE(prob) (prob<=0.0 ? INT_MAX : (int)(-10.0 * log10(prob)))

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

#endif
