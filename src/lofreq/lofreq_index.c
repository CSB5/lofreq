/* -*- c-file-style: "k&r"; indent-tabs-mode: nil; -*- */

/* This is an almost one to one copy of the corresponding bits in samtools */

#include <ctype.h>
#include <assert.h>
#include <stdlib.h>

/* htslib includes */
#include "htslib/faidx.h"
#include "htslib/sam.h"

/* bam_index actually part of API but bam_idxstats not */
int bam_index(int argc, char *argv[]);
int bam_idxstats(int argc, char *argv[]);

/* lofreq includes */
#include "log.h"


#if 1
#define MYNAME "lofreq"
#else
#define MYNAME PACKAGE
#endif


int main_faidx(int argc, char *argv[]) 
{
     char *fa;
     
     fa = argv[2];
     if (fai_build(fa) < 0) {
          return 1;
     }

     return 0;
}

int
main_index(int argc, char *argv[])
{
     char *b = argv[2];
     return sam_index_build(b, 0);
}

int
main_idxstats(int argc, char *argv[])
{
    return bam_idxstats(argc-1, argv+1);
}
