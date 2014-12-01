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
