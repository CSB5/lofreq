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

/* samtools includes */
#include "sam.h"
#include "faidx.h"


/* lofreq includes */
#include "log.h"
#include "utils.h"

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
     int i = -1;
     bam_header_t *header;
     faidx_t *fai;
     char *ref;
     int ref_len = -1;
     bamFile bam_fp;
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

     if (! file_exists(fasta_file)) {
          LOG_FATAL("Fsata file %s does not exist. Exiting...\n", fasta_file);
          return 1;
     }     

     if (0 != strcmp(bam_file, "-")  && ! file_exists(bam_file)) {
          LOG_FATAL("BAM file %s does not exist. Exiting...\n", bam_file);
          return 1;
     }     

     bam_fp = strcmp(bam_file, "-") == 0 ? bam_dopen(fileno(stdin), "r") : bam_open(bam_file, "r");
     header = bam_header_read(bam_fp);
     if (!header) {
          LOG_FATAL("Failed to read BAM header from %s\n", bam_file);
          exit(1);
     }
     
     fai = fai_load(fasta_file);
     
     for (i=0; i < header->n_targets; i++) {
          LOG_DEBUG("BAM header target %d of %d: name=%s len=%d\n", 
                    i+1, header->n_targets,
                    header->target_name[i],
                    header->target_len[i]);
          
          ref = faidx_fetch_seq(fai, header->target_name[i], 
                                0, 0x7fffffff, &ref_len);
          if (NULL == ref) {
               LOG_FATAL("Failed to fetch sequence %s from fasta file\n", header->target_name[i]);
               return -1;
          }
          if (header->target_len[i] != ref_len) {
               LOG_FATAL("Sequence length mismatch for sequence %s (%dbp in fasta; %dbp in bam)\n", 
                         header->target_name[i], header->target_len[i], ref_len);
               return -1;
          }
          free(ref);
     }
     
     fai_destroy(fai);
     bam_header_destroy(header);
     bam_close(bam_fp);

     printf("OK\n");
     return 0;
}
