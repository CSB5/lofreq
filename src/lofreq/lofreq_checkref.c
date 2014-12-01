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


/* This is an almost one to one copy of the corresponding bits in
 * samtools' bam_index.c */

#include <ctype.h>
#include <assert.h>


/* lofreq includes */
#include "log.h"
#include "utils.h"
#include "samutils.h"

#define MYNAME "lofreq checkref"

static void
usage()
{
     fprintf(stderr,
             "\n%s: Check whether given BAM file was created with given reference\n\n", MYNAME);
     fprintf(stderr,"Usage: %s ref.fa in.bam\n\n", MYNAME);
}


int main_checkref(int argc, char *argv[])
{
     char *bam_file;
     char *fasta_file;
     
     if (argc != 4) {
         usage();
         return 1;
     }

     /* get bam file argument
      */
    fasta_file = argv[2];
    bam_file = argv[3];

    if (checkref(fasta_file, bam_file)) {
         printf("Failed\n");
         return 1;
    } else {
         printf("OK\n");
         return 0;
    }
}
