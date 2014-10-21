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

#include <unistd.h>
#include <assert.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <float.h>
#include "faidx.h"
#include "sam.h"
#include "kstring.h"
#include "kaln.h"
#include "kprobaln_ext.h"
#include "defaults.h"

#define USE_EQUAL 1
#define DROP_TAG  2
#define BIN_QUAL  4
#define UPDATE_NM 8
#define UPDATE_MD 16
#define HASH_QNM  32


#define MYNAME "lofreq alnqual"		


#define set_u(u, b, i, k) { int x=(i)-(b); x=x>0?x:0; (u)=((k)-x+1)*3; }
#define prob_to_sangerq(p) (p < 0.0 + DBL_EPSILON ? 126+1 : ((int)(-10 * log10(p))+33))
#define encode_q(q) (uint8_t)(q < 33 ? '!' : (q > 126 ? '~' : q))

/*const char bam_nt16_nt4_table[] = { 4, 0, 1, 4, 2, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4 };*/

int bam_aux_drop_other(bam1_t *b, uint8_t *s);


/* this is lofreq's target function which was heavily modified to accomodate our needs:
 * 1. compute indel alignment qualities on top of base alignment qualities
 * 2. keep base alignment qualities separates, i.e. don't mix with base-qualities
 */
static int bam_prob_realn_core(bam1_t *b, const char *ref, int flag)
{
    int k, i, bw, x, y, z, yb, ye, xb, xe, extend_baq = flag>>1&1;
    /* unused: int apply_baq = flag&1 */
    int redo = 1;  /* = flag&4; */
	uint32_t *cigar = bam1_cigar(b);
	bam1_core_t *c = &b->core;
	kpa_ext_par_t conf = kpa_ext_par_lofreq;
	/*uint8_t *bq = 0, *zq = 0, *qual = bam1_qual(b);*/
	uint8_t *qual = bam1_qual(b);

    /* WARNING/FIXME: redo always on so that baq and ai/ad always get recomputed */

	if ((c->flag & BAM_FUNMAP) || b->core.l_qseq == 0) return -1; // do nothing
	// test if BQ or ZQ is present


    if (redo) {/* if BAQ is present and redo is not set do nothing. note: also means AQ will not be computed */
         uint8_t *bq = NULL; /* pointer to precomputed baq values */
         uint8_t *ai = NULL;
         uint8_t *ad = NULL;

         if ((bq = bam_aux_get(b, BAQ_TAG)) != 0 && *bq == 'Z') {
              bam_aux_del(b, bq);
         }
         if ((ai = bam_aux_get(b, AI_TAG)) != 0 && *ai == 'Z') {
              bam_aux_del(b, ai);
         }
         if ((ad =  bam_aux_get(b, AD_TAG)) != 0 && *ad == 'Z') {
              bam_aux_del(b, ad);
         }
    }


    /* block for reuse of BQ/ZQ deleted */


	// find the start and end of the alignment	
	x = c->pos, y = 0, yb = ye = xb = xe = -1;
	for (k = 0; k < c->n_cigar; ++k) {
		int op, l;
		op = cigar[k]&0xf; l = cigar[k]>>4;
		if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF) {
			if (yb < 0) yb = y;
			if (xb < 0) xb = x;
			ye = y + l; xe = x + l;
			x += l; y += l;
		} else if (op == BAM_CSOFT_CLIP || op == BAM_CINS) y += l;
		else if (op == BAM_CDEL) x += l;
		else if (op == BAM_CREF_SKIP) return -1; // do nothing if there is a reference skip
	}
	// set bandwidth and the start and the end
	bw = 7;
	if (abs((xe - xb) - (ye - yb)) > bw)
		bw = abs((xe - xb) - (ye - yb)) + 3;
	conf.bw = bw;
	xb -= yb + bw/2; if (xb < 0) xb = 0;
	xe += c->l_qseq - ye + bw/2;
	if (xe - xb - c->l_qseq > bw)
		xb += (xe - xb - c->l_qseq - bw) / 2, xe -= (xe - xb - c->l_qseq - bw) / 2;
	{ // glocal
		uint8_t *s, *r, *q, *seq = bam1_seq(b), *bq;
		int *state;
		bq = calloc(c->l_qseq + 1, 1);
		memcpy(bq, qual, c->l_qseq);
		s = calloc(c->l_qseq, 1);
		for (i = 0; i < c->l_qseq; ++i) s[i] = bam_nt16_nt4_table[bam1_seqi(seq, i)];
		r = calloc(xe - xb, 1);
		for (i = xb; i < xe; ++i) {
			if (ref[i] == 0) { xe = i; break; }
			r[i-xb] = bam_nt16_nt4_table[bam_nt16_table[(int)ref[i]]];
		}
		state = calloc(c->l_qseq, sizeof(int));
		q = calloc(c->l_qseq, 1);
          
        double **pd = 0;
        pd = calloc(c->l_qseq+1, sizeof(double*));
          
#ifdef DEBUG
        fprintf(stderr, "processing read %s\n", bam1_qname(b));
#endif
        int bw;
        kpa_ext_glocal(r, xe-xb, s, c->l_qseq, qual, &conf, state, q, pd, &bw);

        /***************************************************************
         * AQ MAGIC START
         */

          // count the number of indels and compute posterior probability
          uint8_t *iaq = 0, *daq = 0;
          int n_ins = 0, n_del = 0;
          iaq = calloc(c->l_qseq + 1, 1);
          daq = calloc(c->l_qseq + 1, 1);

          for (k = 0; k < c->l_qseq; k++) {
               iaq[k] = daq[k] = '~';
          }
          iaq[k] = daq[k] = '\0';
          
          for (k = 0, x = c->pos, y = 0, z = 0; k < c->n_cigar; ++k) { 
               int j, op = cigar[k]&0xf, oplen = cigar[k]>>4;
               // this could be merged into the later block
               if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF) {
                    for (j = 0; j < oplen; j++) {
                         x++; // coordinate on reference
                         y++; // coordinate on query
                         z++; // coordinate on query w/o softclip
                    }
               } else if (op == BAM_CDEL) {
                    if (oplen > 16) continue;
                    n_del += 1;
                    char del_seq[64];
                    int rpos = x; 
                    int qpos = y;
                    if (qpos == 0) continue;
                    for (j = 0; j < oplen; j++) {
                         del_seq[j] = toupper(ref[x]);
                         x++;
                    }
                    del_seq[j] = '\0';
                    int ref_i = x;
                    int del_rep = 0;
                    int rep_i = 0;
                    while (ref_i < xe) {
                         if (ref[ref_i] != del_seq[rep_i]) {
                              break;
                         }
                         del_rep += 1;
                         ref_i += 1;
                         rep_i += 1;
                         if (rep_i >= oplen) {
                              rep_i = 0;
                         }
                    }
                    double ap = 0;
                    for (j = 0; j < del_rep+1; j++) {
                         if (qpos+j > c->l_qseq) break;
                         double *pdi = pd[qpos+j];
                         int u;
                         set_u(u, bw, qpos+j, rpos-xb+1+j);
                         ap += pdi[u+2];
#ifdef DEBUG
                         fprintf(stderr, "probability to add is (%d:%d:%d) %lg\n", 
                             qpos+j, rpos-xb+1+j, u, pdi[u+2]);
#endif
                    }
                    ap = 1 - ap;
                    daq[qpos-1] = encode_q(prob_to_sangerq(ap));
#ifdef DEBUG
                    fprintf(stderr, "DEL %s %d %lg %c %s\n",
                         del_seq, del_rep+1, ap, daq[qpos-1], bam1_qname(b));
#endif
               } else if (op == BAM_CINS) {
                    if (oplen > 16) continue;
                    n_ins += 1;
                    char ins_seq[64];
                    int rpos = x;
                    int qpos = y;
                    if (qpos == 0) continue;
                    for (j = 0; j < oplen; j++) {
                         ins_seq[j] = bam_nt16_rev_table[bam1_seqi(bam1_seq(b), y)];
                         y++;
                         z++;
                    }
                    ins_seq[j] = '\0';
                    int ins_rep = 0;
                    int ref_i = x;
                    int rep_i = 0;
                    while (ref_i < xe) {
                         if (ref[ref_i] != ins_seq[rep_i]) {
                              break;
                         }
                         ins_rep += 1;
                         ref_i += 1;
                         rep_i += 1;
                         if (rep_i >= oplen) {
                              rep_i = 0;
                         }
                    }
                    double ap = 0;
                    for (j = 0; j < ins_rep+1; j++) {
                         if (qpos+j+1 > c->l_qseq) break;
                         double *pdi = pd[qpos+j+1]; 
                         int u;
                         set_u(u, bw, qpos+j+1, rpos-xb+j);
                         ap += pdi[u+1];
#ifdef DEBUG
                         fprintf(stderr, "probability to add is (%d:%d:%d) %lg\n", 
                             qpos+j+1, rpos-xb+j, u, pdi[u+1]);
#endif
                    }
                    ap = 1 - ap; // probability of alignment error
                    iaq[qpos-1] = encode_q(prob_to_sangerq(ap));
#ifdef DEBUG
                    fprintf(stderr, "INS %s %d %lg %c %s\n", 
                         ins_seq, ins_rep+1, ap, iaq[qpos-1], bam1_qname(b));
#endif
               } else if (op == BAM_CSOFT_CLIP) {
                    for (j = 0; j < oplen; j++) {
                         y++;
                    }
               }
          }
 
          for (i = 0; i<=c->l_qseq; ++i) free(pd[i]);
          free(pd); 
          
          /*
           * AQ MAGIC END
           ***************************************************************/
		
          if (!extend_baq) { // in this block, bq[] is capped by base quality qual[]
			for (k = 0, x = c->pos, y = 0; k < c->n_cigar; ++k) {
				int op = cigar[k]&0xf, l = cigar[k]>>4;
				if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF) {
					for (i = y; i < y + l; ++i) {
						if ((state[i]&3) != 0 || state[i]>>2 != x - xb + (i - y)) bq[i] = 0;
#ifdef ORIG_BAQ
						else bq[i] = bq[i] < q[i]? bq[i] : q[i];
#else
                        /* keep the actual values and don't cap by base quality */
                        bq[i] = q[i];
#endif
					}
					x += l; y += l;
				} else if (op == BAM_CSOFT_CLIP || op == BAM_CINS) y += l;
				else if (op == BAM_CDEL) x += l;
			}
#ifdef ORIG_BAQ
			for (i = 0; i < c->l_qseq; ++i) bq[i] = qual[i] - bq[i] + 64; // finalize BQ
#endif

		} else { // in this block, bq[] is BAQ that can be larger than qual[] (different from the above!)
			uint8_t *left, *rght;
			left = calloc(c->l_qseq, 1); rght = calloc(c->l_qseq, 1);
			for (k = 0, x = c->pos, y = 0; k < c->n_cigar; ++k) {
				int op = cigar[k]&0xf, l = cigar[k]>>4;
				if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF) {
					for (i = y; i < y + l; ++i)
						bq[i] = ((state[i]&3) != 0 || state[i]>>2 != x - xb + (i - y))? 0 : q[i];
					for (left[y] = bq[y], i = y + 1; i < y + l; ++i)
						left[i] = bq[i] > left[i-1]? bq[i] : left[i-1];
					for (rght[y+l-1] = bq[y+l-1], i = y + l - 2; i >= y; --i)
						rght[i] = bq[i] > rght[i+1]? bq[i] : rght[i+1];
					for (i = y; i < y + l; ++i)
						bq[i] = left[i] < rght[i]? left[i] : rght[i];
					x += l; y += l;
				} else if (op == BAM_CSOFT_CLIP || op == BAM_CINS) y += l;
				else if (op == BAM_CDEL) x += l;
			}
#ifdef ORIG_BAQ
			for (i = 0; i < c->l_qseq; ++i) bq[i] = 64 + (qual[i] <= bq[i]? 0 : qual[i] - bq[i]); // finalize BQ
#endif
			free(left); free(rght);
		}

#ifndef ORIG_BAQ
          /* need to cap to phred max to be able to store it */
          for (i = 0; i < c->l_qseq; ++i) {
               if (bq[i] > 93) {
                    bq[i] = 93;
               }
               bq[i] += 33;
          }
#endif

#ifdef ORIG_BAQ
		if (apply_baq) {
			for (i = 0; i < c->l_qseq; ++i) qual[i] -= bq[i] - 64; // modify qual
			bam_aux_append(b, "ZQ", 'Z', c->l_qseq + 1, bq);
		} else bam_aux_append(b, "BQ", 'Z', c->l_qseq + 1, bq);
#else
        bam_aux_append(b, BAQ_TAG, 'Z', c->l_qseq + 1, bq);
#endif
        free(bq); free(s); free(r); free(q); free(state);

        /*fprintf(stderr, "%s:%s:%d n_ins=%d n_del=%d\n", __FILE__, __FUNCTION__, __LINE__, n_ins, n_del);*/
        if (n_ins) bam_aux_append(b, AI_TAG, 'Z', c->l_qseq+1, iaq);
        if (n_del) bam_aux_append(b, AD_TAG, 'Z', c->l_qseq+1, daq);
		
        free(iaq); free(daq);
	}
	return 0;
}


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
		/* FIXME: ext should be default */
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
			if (is_realn) bam_prob_realn_core(b, ref, baq_flag);
		}
		samwrite(fpout, b);
	}
	bam_destroy1(b);

	free(ref);
	fai_destroy(fai);
	samclose(fp); samclose(fpout);
	return 0;
}
