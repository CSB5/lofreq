/* -*- c-file-style: "k&r"; indent-tabs-mode: nil; -*- */

/* This is an almost one to one copy of the corresponding bits in samtools */

#include <ctype.h>
#include <assert.h>

/* samtools includes */
#include "bam.h"
#include "sam.h"

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
     char *fi; char *fa;
     
     fa = argv[2];
     fi = samfaipath(fa);
     if (! fi) {
          return 1;
     }

     free(fi);
     return 0;
}

int
main_index(int argc, char *argv[])
{
     char *b = argv[2];
     return bam_index_build(b);
}

int
main_idxstats(int argc, char *argv[])
{
    return bam_idxstats(argc-1, argv+1);
}
