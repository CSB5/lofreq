#ifndef MULTTEST_H
#define MULTTEST_H

/* -*- c-file-style: "k&r"; indent-tabs-mode: nil; -*-
 *
 * FIXME Copyright update
 *
 *********************************************************************/

typedef enum
{
     MTC_NONE,
     MTC_BONF,
     MTC_HOLMBONF,
     MTC_FDR
} mtc_type_t;


#define STR(name) # name

static char *mtc_type_str[] = {
    STR(MTC_NONE),
    STR(MTC_BONF),
    STR(MTC_HOLMBONF),
    STR(MTC_FDR),
};


void
bonf_corr(double data[], int size, int num_tests);

void
holm_bonf_corr(double data[], int size, double alpha, int num_tests);

int
fdr(double data[], int size, double alpha, int num_tests, int **irejected);

int
mtc_str_to_type(char *t);

void
mtc_str(char *buf, int mtc_type);

#endif
