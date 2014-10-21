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

#include <stdlib.h>
#include "faidx.h"
#include "sam.h"
#include "bam_md_ext.h"
#include "defaults.h"


#define USE_EQUAL 1
#define DROP_TAG  2
#define BIN_QUAL  4
#define UPDATE_NM 8
#define UPDATE_MD 16
#define HASH_QNM  32

#define MYNAME "lofreq alnqual"		



int main(int argc, char *argv[])
{
	int c, flt_flag, tid = -2, ret, len, is_bam_out, is_sam_in, is_uncompressed, max_nm, is_realn, capQ, baq_flag;
	samfile_t *fp, *fpout = 0;
	faidx_t *fai;
	char *ref = 0, mode_w[8], mode_r[8];
	bam1_t *b;

	flt_flag = UPDATE_NM | UPDATE_MD;
	is_bam_out = is_sam_in = is_uncompressed = is_realn = max_nm = capQ = baq_flag = 0;
	mode_w[0] = mode_r[0] = 0;
	strcpy(mode_r, "r"); strcpy(mode_w, "w");
	
	is_realn = 1;
	baq_flag |= 2; /* ext baq */
	while ((c = getopt(argc, argv, "buSe")) >= 0) {
		switch (c) {
#if 0
		case 'r': is_realn = 1; break;
		case 'e': flt_flag |= USE_EQUAL; break;
		case 'A': baq_flag |= 1; break;
		case 'd': flt_flag |= DROP_TAG; break;
		case 'q': flt_flag |= BIN_QUAL; break;
		case 'h': flt_flag |= HASH_QNM; break;
		case 'N': flt_flag &= ~(UPDATE_MD|UPDATE_NM); break;
		case 'n': max_nm = atoi(optarg); break;
		case 'C': capQ = atoi(optarg); break;
#endif
		case 'b': is_bam_out = 1; break;
		case 'u': is_uncompressed = is_bam_out = 1; break;
		case 'S': is_sam_in = 1; break;
		case 'e': baq_flag &= ~2; break;
		default: fprintf(stderr, "%s unrecognized option '-%c'\n", MYNAME, c); return 1;
		}
	}
	if (!is_sam_in) strcat(mode_r, "b");
	if (is_bam_out) strcat(mode_w, "b");
	else strcat(mode_w, "h");
	if (is_uncompressed) strcat(mode_w, "u");
	if (optind + 1 >= argc) {
#if 0
  		/* donwannah */
		fprintf(stderr, "Usage:   samtools calmd [-eubrS] <aln.bam> <ref.fasta>\n\n");
		fprintf(stderr, "Options: -e       change identical bases to '='\n");		
		fprintf(stderr, "         -A       modify the quality string\n");
		/* default */
		fprintf(stderr, "         -r       compute the BQ tag (without -A) or cap baseQ by BAQ (with -A)\n");
#else
		fprintf(stderr, "%s: add base- and indel-alignment qualities (BAQ, IDAQ) to BAM file\n\n", MYNAME);
		fprintf(stderr, "Usage:   %s [options] <aln.bam> <ref.fasta>\n", MYNAME);
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "         -u       uncompressed BAM output (for piping)\n");
		fprintf(stderr, "         -b       compressed BAM output\n");
		fprintf(stderr, "         -S       the input is SAM with header\n");
		fprintf(stderr, "         -e       use default instead of extended BAQ (the latter gives better sensitivity but lower specificity)\n\n");		
		fprintf(stderr, "- Output will be written to stdout.\n");				
		fprintf(stderr, "- Only reads containing indels will contain indel-alignment qualities (tags: %s and %s).\n", AI_TAG, AD_TAG);
		fprintf(stderr, "- Do not change the alignmnent after running this, i.e. use this as last postprocessing step!\n");
		fprintf(stderr, "- This program is based on samtools. BAQ was introduced by Heng Li PMID:21320865\n\n");
#endif
		return 1;
	}
	fp = samopen(argv[optind], mode_r, 0);
	if (fp == 0) return 1;
	if (is_sam_in && (fp->header == 0 || fp->header->n_targets == 0)) {
         fprintf(stderr, "%s: input SAM does not have header. Abort!\n", MYNAME);
		return 1;
	}
	fpout = samopen("-", mode_w, fp->header);
	fai = fai_load(argv[optind+1]);

	b = bam_init1();
	while ((ret = samread(fp, b)) >= 0) {
		if (b->core.tid >= 0) {
			if (tid != b->core.tid) {
				free(ref);
				ref = fai_fetch(fai, fp->header->target_name[b->core.tid], &len);
				tid = b->core.tid;
				if (ref == 0)
					fprintf(stderr, "%s fail to find sequence '%s' in the reference.\n",
							MYNAME, fp->header->target_name[tid]);
			}
			if (is_realn) bam_prob_realn_core_ext(b, ref, baq_flag);
		}
		samwrite(fpout, b);
	}
	bam_destroy1(b);

	free(ref);
	fai_destroy(fai);
	samclose(fp); samclose(fpout);
	return 0;
}
