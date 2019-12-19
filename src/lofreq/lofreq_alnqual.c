/* -*- c-file-style: "k&r"; indent-tabs-mode: nil; -*- */
/*
  This is part of LoFreq Star and largely based on samtools' bam_md.c
  (0.1.19) which was originally published under the MIT License.
  
  Copyright (c) 2003-2006, 2008-2010, by Heng Li <lh3lh3@live.co.uk>
  
  Permission is hereby granted, free of charge, to any person obtaining
  a copy of this software and associated documentation files (the
  "Software"), to deal in the Software without restriction, including
  without limitation the rights to use, copy, modify, merge, publish,
  distribute, sublicense, and/or sell copies of the Software, and to
  permit persons to whom the Software is furnished to do so, subject to
  the following conditions:
  
  The above copyright notice and this permission notice shall be
  included in all copies or substantial portions of the Software.
  
  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
  EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
  MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
  NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
  BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
  ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
  CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
  SOFTWARE.
*/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include "htslib/faidx.h"
#include "htslib/sam.h"

#include "utils.h"
#include "bam_md_ext.h"
#include "defaults.h"

#define USE_EQUAL 1
#define DROP_TAG  2
#define BIN_QUAL  4
#define UPDATE_NM 8
#define UPDATE_MD 16
#define HASH_QNM  32

#define MYNAME "lofreq alnqual"		


static void usage()
{
     fprintf(stderr, "%s: add base- and indel-alignment qualities (BAQ, IDAQ) to BAM file\n\n", MYNAME);
     fprintf(stderr, "Usage:   %s [options] <aln.bam> <ref.fasta>\n", MYNAME);
     fprintf(stderr, "Options:\n");
     fprintf(stderr, "         -b       BAM output (instead of SAM)\n");
     fprintf(stderr, "         -u       Uncompressed BAM output (for piping)\n");
     fprintf(stderr, "         -S       The input is SAM with header\n");
     fprintf(stderr, "         -e       Use default instead of extended BAQ (the latter gives better sensitivity but lower specificity)\n");
     fprintf(stderr, "         -B       Don't compute base alignment qualities\n");
     fprintf(stderr, "         -A       Don't compute indel alignment qualities\n");
     fprintf(stderr, "         -r       Recompute i.e. overwrite existing values\n");
     fprintf(stderr, "- Output BAM will be written to stdout.\n");				
     fprintf(stderr, "- Only reads containing indels will contain indel-alignment qualities (tags: %s and %s).\n", AI_TAG, AD_TAG);
     fprintf(stderr, "- Do not change the alignmnent after running this, i.e. use this as last postprocessing step!\n");
     fprintf(stderr, "- This program is based on samtools. BAQ was introduced by Heng Li PMID:21320865\n\n");

}


int main_alnqual(int argc, char *argv[])
{
     int c, tid = -2, ret, len, is_bam_out, is_sam_in, is_uncompressed;
     samFile *fp, *fpout = 0;
     faidx_t *fai;
     char *ref = 0, mode_w[8], mode_r[8];
     bam1_t *b;
     int baq_flag = 1;
     int ext_baq = 1;
     int idaq_flag = 1;
     int redo = 0;

     is_bam_out = is_sam_in = is_uncompressed = 0;
     mode_w[0] = mode_r[0] = 0;
     strcpy(mode_r, "r"); strcpy(mode_w, "w");
	
     while ((c = getopt(argc, argv, "buSeBAr")) >= 0) {
          switch (c) {
          case 'b': is_bam_out = 1; break;
          case 'u': is_uncompressed = is_bam_out = 1; break;
          case 'S': is_sam_in = 1; break;
          case 'e': ext_baq = 0; break;
          case 'B': baq_flag = 0; break;
          case 'A': idaq_flag = 0; break;
          case 'r': redo = 1; break;
          case '?': 
               fprintf(stderr, "FATAL: unrecognized arguments found. Exiting...\n");
               return 1;
          default: 
               break;
          }
     }
     if (optind + 1 >= argc) {
          usage();
          return 1;
     }

     if (!is_sam_in) strcat(mode_r, "b");
     if (is_bam_out) {
          strcat(mode_w, "b");
     } else{
          strcat(mode_w, "h");
     }
     if (is_uncompressed) strcat(mode_w, "u");
     
     if (redo) {
          if (baq_flag) {
               baq_flag = 2;
          }
          if (idaq_flag) {
               idaq_flag = 2;
          }
     }

     if (! baq_flag && ! idaq_flag) {
          fprintf(stderr, "FATAL: %s: Nothing to do: BAQ and IDAQ off\n", MYNAME); 
          return 1;
     }

     fp = sam_open(argv[optind], mode_r);
     if (fp == 0) return 1;
     bam_hdr_t *header = sam_hdr_read(fp);
     if (header == 0) {
          fprintf(stderr, "FATAL: %s: input SAM does not have header\n", MYNAME);
          return 1;
     }
     fpout = sam_open("-", mode_w);
     if (sam_hdr_write(fpout, header) < 0) {
          fprintf(stderr, "FATAL: %s: failed to copy SAM header to output\n", MYNAME);
          return 1;
     }

     fai = fai_load(argv[optind+1]);
     if (! fai) {
          fprintf(stderr, "FATAL: %s: failed to load fai index\n", MYNAME);
          return 1;
     }

     b = bam_init1();
     while ((ret = sam_read1(fp, header, b)) >= 0) {
          if (b->core.tid >= 0) {
               if (tid != b->core.tid) {
                    free(ref);
                    ref = fai_fetch(fai, header->target_name[b->core.tid], &len);
                    strtoupper(ref);/* safeguard */
                    tid = b->core.tid;
                    if (ref == 0) {
                         fprintf(stderr, "FATAL: %s failed to find sequence '%s' in the reference.\n",
                                   MYNAME, header->target_name[tid]);
                         return 1;
                    }
               }
               
               bam_prob_realn_core_ext(b, ref, baq_flag, ext_baq, idaq_flag);
          }
          if (sam_write1(fpout, header, b) < 0) {
                         fprintf(stderr, "FATAL: %s failed to write record to output\n",
                                   MYNAME);
                         return 1;
          }
     }
     bam_destroy1(b);
     
     free(ref);
     fai_destroy(fai);
     bam_hdr_destroy(header);
     sam_close(fp);
     sam_close(fpout);
     return 0;
}
