#ifndef MULTTEST_H
#define MULTTEST_H

/* -*- c-file-style: "k&r"; indent-tabs-mode: nil; -*- */

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
bonf_corr(double data[], long int size, long int num_tests);

void
holm_bonf_corr(double data[], long int size, double alpha, long int num_tests);

long int
fdr(double data[], long int size, double alpha, long int num_tests, long int **irejected);

int
mtc_str_to_type(char *t);

void
mtc_str(char *buf, int mtc_type);

#endif
