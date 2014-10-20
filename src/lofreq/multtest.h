/* -*- c-file-style: "k&r"; indent-tabs-mode: nil; -*- */
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
#ifndef MULTTEST_H
#define MULTTEST_H


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
