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
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <float.h>

#include "htslib/sam.h"
#include "htslib/faidx.h"
#include "htslib/kstring.h"

#include "kprobaln_ext.h"
#include "samutils.h"
#include "defaults.h"
#include "bam_md_ext.h"

#ifdef PACBIO_REALN
static int pacbio_msg_printed = 0;
#endif



void idaq(bam1_t *b, const char *ref, double **pd, int xe, int xb, int bw);

#define set_u(u, b, i, k) { int x=(i)-(b); x=x>0?x:0; (u)=((k)-x+1)*3; }
#define prob_to_sangerq(p) (p < 0.0 + DBL_EPSILON ? 126+1 : ((int)(-10 * log10(p))+33))
#define encode_q(q) (uint8_t)(q < 33 ? '!' : (q > 126 ? '~' : q))


/* fw and bck matrices in kprob have alloc limit
   bw2 = bw * 2 + 1
   alloc(bw2 * 3 + 6)
   and in addition the original BAQ checks whether if (u < 3 || u >= bw2*3+3) and continues if so
*/
int u_within_limits(int u, int bw) {
     int bw2 = bw * 2 + 1;
     if (u<3 || u >= bw2*3+3) {
          return 0;
     } else {
          return 1;
     }
}     

void idaq(bam1_t *b, const char *ref, double **pd, int xe, int xb, int bw)
{
	uint32_t *cigar = bam_get_cigar(b);
	bam1_core_t *c = &b->core;
    // count the number of indels and compute posterior probability
    uint8_t *iaq = 0, *daq = 0;
    int n_ins = 0, n_del = 0;
    int k, x, y, z;

#if 0
    fprintf(stderr, "Running idaq on %s with cigar %s\n", bam_get_qname(b), cigar_str_from_bam(b));
#endif

    iaq = calloc(c->l_qseq + 1, 1);
    daq = calloc(c->l_qseq + 1, 1);
    
    /* init to highest possible value */
    for (k = 0; k < c->l_qseq; k++) {
         iaq[k] = daq[k] = '~';
    }
    iaq[k] = daq[k] = '\0';
    
    /* equivalent indels may occur in repetitive regions. In such
     * cases, we estimate the alignment probability of an indel event
     * as the sum of the alignment probability of all equivalent indel
     * events. see del_rep and ins_rep handling below 
     */
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
              char *del_seq;
              int rpos = x; 
              int qpos = y;
              int ref_i;
              int del_rep = 0;/* if in repetetive region */
              int rep_i = 0;
              double ap = 0;

              if (qpos == 0) continue;
              if (oplen > 16) continue; /*FIXME why */
              n_del += 1;
              del_seq = malloc((oplen+1)*sizeof(char));
              for (j = 0; j < oplen; j++) {
                   del_seq[j] = ref[x];
                   x++;
              }
              del_seq[j] = '\0';
              ref_i = x;
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
              for (j = 0; j < del_rep+1; j++) {
                   if (qpos+j > c->l_qseq) break;
                   double *pdi = pd[qpos+j];
                   int u;

                   set_u(u, bw, qpos+j, rpos-xb+1+j);
                   /* FIXME happens for long reads, i.e. pacbio. why? see corresponding bit for ins_rep 
                    */
                   if (! u_within_limits(u, bw)) {
#if 0
                        fprintf(stderr, "WARNING u of %d not within limits for %s\n", u, bam_get_qname(b));
#endif
                        continue;
                   }
                   ap += pdi[u+2];
#if 0
                   fprintf(stderr, "probability to add comes from pd[%d+%d + %d+%d = %d]. qseq+1 is %d\n", 
                           qpos,j,u,2, qpos+j+u+2, c->l_qseq+1);
                   fprintf(stderr, "probability to add is (%d:%d:%d) %lg\n", 
                           qpos+j, rpos-xb+1+j, u, pdi[u+2]);
                   fflush(stderr);

#endif
              }
              ap = 1 - ap;
              daq[qpos-1] = encode_q(prob_to_sangerq(ap));
              /*fprintf(stderr, "DAQ %d: %c %g\n", qpos-1, daq[qpos-1], ap);*/
              free(del_seq);
#ifdef DEBUG
              fprintf(stderr, "DEL %s %d %lg %c %s\n",
                      del_seq, del_rep+1, ap, daq[qpos-1], bam_get_qname(b));
#endif
         } else if (op == BAM_CINS) {
              char *ins_seq;
              int rpos = x;
              int qpos = y;
              int ins_rep = 0; /* if in repetetive region */
              int ref_i = x;
              int rep_i = 0;
              double ap = 0;

              if (oplen > 16) continue; /*FIXME why */
              n_ins += 1;
              if (qpos == 0) continue;
              ins_seq = malloc((oplen+1)*sizeof(char));
              for (j = 0; j < oplen; j++) {
                   ins_seq[j] = seq_nt16_str[bam_seqi(bam_get_seq(b), y)];
                   y++;
                   z++;
              }
              ins_seq[j] = '\0';
              ref_i = x;
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
              for (j = 0; j < ins_rep+1; j++) {
                   if (qpos+j+1 > c->l_qseq) break;
                   double *pdi = pd[qpos+j+1]; 
                   int u;

                   set_u(u, bw, qpos+j+1, rpos-xb+j);
                   /* FIXME happens for long reads, i.e. pacbio. why? see corresponding bit for del_rep 
                    */
                   if (! u_within_limits(u, bw)) {
#if 0
                        fprintf(stderr, "WARNING u of %d not within limits for %s\n", u, bam_get_qname(b));
#endif
                        continue;
                   }
                   ap += pdi[u+1];
#if 0
                   fprintf(stderr, "probability to add comes from pd[%d+%d+%d + %d+%d = %d]. qseq+1 is %d\n", 
                           qpos,j,1,u,1, qpos+j+1+u+1, c->l_qseq+1);
                   fprintf(stderr, "probability to add is (%d:%d:%d) %lg\n", 
                           qpos+j+1, rpos-xb+j, u, pdi[u+1]);
                   fflush(stderr);
#endif
              }
              ap = 1 - ap; // probability of alignment error
              iaq[qpos-1] = encode_q(prob_to_sangerq(ap));
              /*fprintf(stderr, "IAQ %d: %c %g\n", qpos-1, iaq[qpos-1], ap);*/
              free(ins_seq);
#ifdef DEBUG
              fprintf(stderr, "INS %s %d %lg %c %s\n", 
                      ins_seq, ins_rep+1, ap, iaq[qpos-1], bam_get_qname(b));
#endif
         } else if (op == BAM_CSOFT_CLIP) {
              for (j = 0; j < oplen; j++) {
                   y++;
              }
         }
    }
    
    /*fprintf(stderr, "%s:%s:%d n_ins=%d n_del=%d\n", __FILE__, __FUNCTION__, __LINE__, n_ins, n_del);*/
    if (n_ins) {
         bam_aux_append(b, AI_TAG, 'Z', c->l_qseq+1, iaq);
    }
    if (n_del)  {
         bam_aux_append(b, AD_TAG, 'Z', c->l_qseq+1, daq);
    }

    free(iaq); free(daq);
}



/* this is lofreq's target function which was heavily modified to accomodate our needs:
 * 1. compute indel alignment qualities on top of base alignment qualities
 * 2. keep base alignment qualities separates, i.e. don't mix with base-qualities
 *
 * baq_flag: 0 off, 1 on, 2 redo
 * aq_flag: 0 off, 1 on, 2 redo
 */
int bam_prob_realn_core_ext(bam1_t *b, const char *ref, 
                            int baq_flag, int baq_extended,
                            int idaq_flag)
{
/*#define ORIG_BAQ 1*/
     int k, i, bw, x, y, yb, ye, xb, xe;
     uint32_t *cigar = bam_get_cigar(b);
     bam1_core_t *c = &b->core;
#ifdef PACBIO_REALN
     kpa_ext_par_t conf = kpa_ext_par_lofreq_pacbio;
     if (! pacbio_msg_printed) {
          fprintf(stderr, "WARN(%s|%s): Using pacbio viterbi params\n", __FILE__, __FUNCTION__);
          pacbio_msg_printed = 1;
     }
#else
     kpa_ext_par_t conf = kpa_ext_par_lofreq_illumina;
#endif
     /*uint8_t *bq = 0, *zq = 0, *qual = bam_get_qual(b);*/
     uint8_t *qual = bam_get_qual(b);
     uint8_t *prec_ai, *prec_ad, *prec_baq;
     int has_ins = 0, has_del = 0;
     double **pd = 0;

     /* nothing to do ? */
     if (! baq_flag && ! idaq_flag) {
          return 0;
     }

     /*fprintf(stderr, "FIXME baq_flag=%d idaq_flag=%d\n", baq_flag, idaq_flag);*/

     /* no alignment? */
     if ((c->flag & BAM_FUNMAP) || b->core.l_qseq == 0) {
          return 0;
     }
     
     /* get existing tags. delete if existing and redo is on
      */
     if ((prec_baq = bam_aux_get(b, BAQ_TAG)) != 0 && *prec_baq == 'Z') {
          if (baq_flag==2) {
               bam_aux_del(b, prec_baq);
               prec_baq = NULL;
          }
     }
     if ((prec_ai = bam_aux_get(b, AI_TAG)) != 0 && *prec_ai == 'Z') {
          if (idaq_flag==2) {
               bam_aux_del(b, prec_ai);
               prec_ai = NULL;
          }
     }
     if ((prec_ad = bam_aux_get(b, AD_TAG)) != 0 && *prec_ad == 'Z') {
          if (idaq_flag==2) {
               bam_aux_del(b, prec_ad);
               prec_ad = NULL;
          }
     }

	/* find the start and end of the alignment */
	x = c->pos, y = 0, yb = ye = xb = xe = -1;
	for (k = 0; k < c->n_cigar; ++k) {
		int op, l;
		op = cigar[k]&0xf; l = cigar[k]>>4;
		if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF) {
			if (yb < 0) yb = y;
			if (xb < 0) xb = x;
			ye = y + l; xe = x + l;
			x += l; y += l;
		} else if (op == BAM_CSOFT_CLIP || op == BAM_CINS) {
             y += l;
             if (op == BAM_CINS) {
                  has_ins = 1;
             }
		} else if (op == BAM_CDEL) {
             has_del = 1;
             x += l;
        }
		else if (op == BAM_CREF_SKIP) {
#if 0
             return 0; /* do nothing if there is a reference skip */
#else
             /* returning would mean give up and compute no BAQ. 
                behaviour now modelled after calc_read_alnerrprof(),
                where CDEL and CREF_SKIP behave the same */
             x += l; 
#endif
        }
	}

#if 0
    fprintf(stderr, "%s with cigar %s: baq_flag=%d prec_baq=%p has_del=%d prec_ad=%p has_ins=%d prec_ai=%p, idaq_flag=%d\n", 
            bam_get_qname(b), cigar_str_from_bam(b),  baq_flag, prec_baq, has_del, prec_ad, has_ins, prec_ai, idaq_flag);
#endif
    /* don't do anything if everything's there already */
    if (baq_flag==0 || prec_baq) {
         int skip = 1;
         if (has_del && ! prec_ad) {
              skip = 0;
         }
         if (has_ins && ! prec_ai) {
              skip = 0;
         }
         if (skip) {
#if 0
              fprintf(stderr, "Reusing all alignment quality values for read %s!\n", bam_get_qname(b));
#endif
              return 0;
         }
    }

    if (has_ins || has_del) {
         pd = calloc(c->l_qseq+1, sizeof(double*));
    }

    /* either need to compute BAQ or IDAQ 
     */

	/* set bandwidth and the start and the end */
	bw = 7;
	if (abs((xe - xb) - (ye - yb)) > bw)
		bw = abs((xe - xb) - (ye - yb)) + 3;
	conf.bw = bw;
	xb -= yb + bw/2; if (xb < 0) xb = 0;
	xe += c->l_qseq - ye + bw/2;
	if (xe - xb - c->l_qseq > bw)
		xb += (xe - xb - c->l_qseq - bw) / 2, xe -= (xe - xb - c->l_qseq - bw) / 2;


	{ /* glocal */
		uint8_t *s, *r, *q, *seq = bam_get_seq(b), *bq;
		int *state;
        int bw;

		bq = calloc(c->l_qseq + 1, 1);
		memcpy(bq, qual, c->l_qseq);
		s = calloc(c->l_qseq, 1);
		for (i = 0; i < c->l_qseq; ++i) s[i] = seq_nt16_int[bam_seqi(seq, i)];
		r = calloc(xe - xb, 1);
		for (i = xb; i < xe; ++i) {
			if (ref[i] == 0) { xe = i; break; }
			r[i-xb] = seq_nt16_int[seq_nt16_table[(int)ref[i]]];
		}
		state = calloc(c->l_qseq, sizeof(int));
		q = calloc(c->l_qseq, 1);
          
          
#ifdef DEBUG
        fprintf(stderr, "processing read %s\n", bam_get_qname(b));
#endif
        kpa_ext_glocal(r, xe-xb, s, c->l_qseq, qual, &conf, state, q, pd, &bw);

        if (baq_flag && ! prec_baq) {
             if (! baq_extended) { // in this block, bq[] is capped by base quality qual[]
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
                  if (bq[i] > SANGER_PHRED_MAX) {
                       bq[i] = SANGER_PHRED_MAX;
                  }
                  bq[i] += 33;
             }
#endif
             
/*#undef ORIG_BAQ*/
#ifdef ORIG_BAQ
             if (apply_baq) {
                  for (i = 0; i < c->l_qseq; ++i) qual[i] -= bq[i] - 64; // modify qual
                  bam_aux_append(b, "ZQ", 'Z', c->l_qseq + 1, bq);
             } else bam_aux_append(b, "BQ", 'Z', c->l_qseq + 1, bq);
#else
             bam_aux_append(b, BAQ_TAG, 'Z', c->l_qseq + 1, bq);
#endif
        }
        /* no baq */
        
        
        if (idaq_flag && pd) {/* pd served as previous check to see if ai or ad actually need to be computed */
             idaq(b, ref, pd, xe, xb, bw);
        }
        
        if (pd) {
             for (i = 0; i<=c->l_qseq; ++i) free(pd[i]);
             free(pd); 
        }
        free(bq); free(s); free(r); free(q); free(state);
	}

	return 0;
}

