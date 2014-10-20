/* -*- c-file-style: "k&r"; indent-tabs-mode: nil; -*- */

/* This is an almost one to one copy of the corresponding bits in
 * samtools' bam_index.c */

#include <ctype.h>
#include <assert.h>

/* samtools includes */
#include "bam.h"
int bam_index(int argc, char *argv[]);
int bam_idxstats(int argc, char *argv[]);

/* lofreq includes */
#include "log.h"


#if 1
#define MYNAME "lofreq"
#else
#define MYNAME PACKAGE
#endif


int
main_index(int argc, char *argv[])
{
    return bam_index(argc-1, argv+1);
}

int
main_idxstats(int argc, char *argv[])
{
    return bam_idxstats(argc-1, argv+1);
}
