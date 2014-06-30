/* -*- c-file-style: "k&r"; indent-tabs-mode: nil; -*- */

/*********************************************************************
 *
 * FIXME update license
 *
 *********************************************************************/

/* This is an almost one to one copy of the corresponding bits in
 * samtools' bam_index.c */

#include <ctype.h>
#include <assert.h>


/* lofreq includes */
#include "log.h"
#include "utils.h"
#include "samutils.h"

#define MYNAME "lofreq checkref"

void
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
