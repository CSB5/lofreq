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


#include <ctype.h>
#include <stdio.h>
#include <getopt.h>
#include <stdlib.h>

#include "htslib/faidx.h"
#include "sam.h"
#include "viterbi.h"
#include "log.h"
#include "lofreq_viterbi.h"
#include "utils.h"

#define SANGERQUAL_TO_PHRED(c) ((int)(c)-33)

/* FIXME: implement auto clipping of Q2 tails */

#define RWIN 10

typedef struct {
     samfile_t *in;
     bamFile out;
     faidx_t *fai;
     uint32_t tid;
     char *ref;
     int reflen;
} tmpstruct_t;

static void replace_cigar(bam1_t *b, int n, uint32_t *cigar)
{
    if (n != b->core.n_cigar) {
        int o = b->core.l_qname + b->core.n_cigar * 4;
        if (b->data_len + (n - b->core.n_cigar) * 4 > b->m_data) {
            b->m_data = b->data_len + (n - b->core.n_cigar) * 4;
            kroundup32(b->m_data);
            b->data = (uint8_t*)realloc(b->data, b->m_data);
        }
        memmove(b->data + b->core.l_qname + n * 4, b->data + o, b->data_len - o);
        memcpy(b->data + b->core.l_qname, cigar, n * 4);
        b->data_len += (n - b->core.n_cigar) * 4;
        b->core.n_cigar = n;
    } else memcpy(b->data + b->core.l_qname, cigar, n * 4);
}


/* function checks if alignment is made of all Q2s
 * if not, returns remaining values so that median 
 *  can be calculated
 */
int check_Q2(char *bqual, int *num){
    int is_all_Q2 = 1;
    int i, pom = 0;
    int l = strlen(bqual);
    *num = 0;
    for (i=0; i<l; i++){
            if (SANGERQUAL_TO_PHRED(bqual[i]) != 2){
                    pom++;
                    is_all_Q2 = 0;
            }
    }       
    *num = pom;
    return is_all_Q2;
}

void remain(char *bqual, int *remaining){
     int pom = 0;
     int i, q;
     int l = strlen(bqual);
     for (i=0; i<l; i++){
          q = SANGERQUAL_TO_PHRED(bqual[i]);
          if (q != 2){
               remaining[pom] = q;
               pom++;
          }
     }   
}

static int fetch_func(bam1_t *b, void *data, int del_flag, int q2def, int reclip)
{
     /* see
      https://github.com/lh3/bwa/blob/426e54740ca2b9b08e013f28560d01a570a0ab15/ksw.c
      for optimizations and speedups
     */
     tmpstruct_t *tmp = (tmpstruct_t*)data;
     bam1_core_t *c = &b->core;
     uint8_t *seq = bam1_seq(b);
     uint32_t *cigar = bam1_cigar(b);
     int reflen;
    
     if (del_flag) {
          uint8_t *old_nm;
          uint8_t *old_mc;
          uint8_t *old_md;
          uint8_t *old_as;

          /* once you bam_aux_del b will change and all pointers to it, so don't use bam_aux_get again too early */
          
          old_nm = bam_aux_get(b, "NM");          
          if (old_nm) {          
               bam_aux_del(b, old_nm);
          }

          old_mc = bam_aux_get(b, "MC");          
          if (old_mc) {          
                bam_aux_del(b, old_mc);          
          }

          old_md = bam_aux_get(b, "MD");          
          if (old_md) {          
               bam_aux_del(b, old_md);
          }
          
          old_as = bam_aux_get(b, "AS");                    
          if (old_as) {
               bam_aux_del(b, old_as);
          }
     }

     if (c->flag & BAM_FUNMAP) {
          bam_write1(tmp->out, b);
          return 0;
     }

     /* fetch reference sequence if incorrect tid */
     if (tmp->tid != c->tid) {
          if (tmp->ref) free(tmp->ref);
          if ((tmp->ref = 
               fai_fetch(tmp->fai, tmp->in->header->target_name[c->tid], &reflen)) == 0) {
               fprintf(stderr, "failed to find reference sequence %s\n", 
                                tmp->in->header->target_name[c->tid]);
          }
          strtoupper(tmp->ref);/* safeguard */
          tmp->tid = c->tid;
          tmp->reflen = reflen;
     }
     int i;

     // remove soft clipped bases
     char query[c->l_qseq+1];
     char bqual[c->l_qseq+1];

     int x = c->pos; // coordinate on reference
     int y = 0; // coordinate on query
     int z = 0; // coordinate on query w/o softclip

     int indels = 0;

     // parse cigar string
     for (i = 0; i < c->n_cigar; ++i) {
          int j, oplen = cigar[i] >> 4, op = cigar[i]&0xf;
          if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF) {
               for (j = 0; j < oplen; j++) {
                    query[z] = bam_nt16_rev_table[bam1_seqi(seq, y)];
                    bqual[z] = (char)bam1_qual(b)[y]+33;
                    x++;
                    y++;
                    z++;
               }
          } else if (op == BAM_CHARD_CLIP) {
               /* in theory we should do nothing here but hard clipping info gets lost here FIXME
                */               
               bam_write1(tmp->out, b);
               return 1;
          } else if (op == BAM_CDEL) {
               x += oplen;
               indels += 1;
          } else if (op == BAM_CINS) {
               for (j = 0; j < oplen; j++) {
                    query[z] = bam_nt16_rev_table[bam1_seqi(seq, y)];
                    bqual[z] = (char)bam1_qual(b)[y]+33;
                    y++;
                    z++;
               }
               indels += 1;
          } else if (op == BAM_CSOFT_CLIP) {
               for (j = 0; j < oplen; j++) {
                    y++;
               }
          } else {
               LOG_WARN("Unknown cigar op %d. Not touching read %s\n", op, bam1_qname(b));
               bam_write1(tmp->out, b);
               return 1;
          }
     }
     query[z] = bqual[z] = '\0';

     if (indels == 0) {
          bam_write1(tmp->out, b);
          return 0;
     }
    int len_remaining = 0;
    if (check_Q2(bqual, &len_remaining)) {
		if (reclip){
			// check if first op or last op is I and replace with S
			 int curr_oplen_check = cigar[0] >> 4;
			 int curr_op_check = cigar[0]&0xf;
			 if (curr_op_check == BAM_CINS){
				curr_op_check = BAM_CSOFT_CLIP;
				cigar[0] = curr_oplen_check <<4 | curr_op_check;
			}
			 curr_oplen_check = cigar[c->n_cigar-1] >> 4;
			 curr_op_check = cigar[c->n_cigar-1]&0xf;

			 if (curr_op_check == BAM_CINS){
				curr_op_check = BAM_CSOFT_CLIP;
				cigar[c->n_cigar-1] = curr_oplen_check <<4 | curr_op_check;
			}
			
			replace_cigar(b,c->n_cigar,cigar);
		}
        bam_write1(tmp->out, b);
        return 0;
    }
    int remaining[len_remaining+1];
    remain(bqual, remaining);
    remaining[len_remaining] = '\0';
    if (q2def < 0) {
        q2def = int_median(remaining, len_remaining);
    }
    
     /* get reference with RWIN padding */
     char ref[c->l_qseq+1+indels+RWIN*2];
     int lower = c->pos - RWIN;
     lower = lower < 0? 0: lower;
     int upper = x + RWIN;
     upper = upper > tmp->reflen? tmp->reflen: upper;
     for (z = 0, i = lower; i < upper; z++, i++) {
          ref[z] = tmp->ref[i];
     }
     ref[z] = '\0';

     /* run viterbi */
     char *aln = malloc(sizeof(char)*(2*(c->l_qseq)));
     int shift = viterbi(ref, query, bqual, aln, q2def);

     /* convert to cigar */
     uint32_t *realn_cigar = 0;
     int realn_n_cigar = 0;
     
     /* check if soft-clipped in the front */
     int curr_oplen = cigar[0] >> 4; 
     int curr_op = cigar[0]&0xf;
     if (curr_op == BAM_CSOFT_CLIP) {
          realn_cigar = realloc(realn_cigar, (realn_n_cigar+1)*sizeof(uint32_t));
          realn_cigar[realn_n_cigar] = curr_oplen<<4 | curr_op;
          realn_n_cigar += 1;
     }
     
     /* get cigar of the realigned query */
     curr_op = aln[0] == 'M' ? 0 : (aln[0] == 'I'? 1 : 2);
     curr_oplen = 1;
     for (i = 1; i < strlen(aln); i++) {
          int this_op = aln[i] == 'M' ? 0 : (aln[i] == 'I' ? 1 : 2);
          if (this_op != curr_op) {
               realn_cigar = realloc(realn_cigar, (realn_n_cigar+1)*sizeof(uint32_t));
               realn_cigar[realn_n_cigar] = curr_oplen<<4 | curr_op;
               realn_n_cigar += 1;
               curr_op = this_op;
               curr_oplen = 1;
          } else {
               curr_oplen += 1;
          }
     }
     realn_cigar = realloc(realn_cigar, (realn_n_cigar+1)*sizeof(uint32_t));
     realn_cigar[realn_n_cigar] = curr_oplen<<4 | curr_op;
     realn_n_cigar += 1; 
    
     /* check if soft-clipped in the back */
     curr_oplen = cigar[c->n_cigar-1] >> 4; 
     curr_op = cigar[c->n_cigar-1]&0xf;
     if (curr_op == BAM_CSOFT_CLIP) {
          realn_cigar = realloc(realn_cigar, (realn_n_cigar+1)*sizeof(uint32_t));
          realn_cigar[realn_n_cigar] = curr_oplen<<4 | curr_op;
          realn_n_cigar += 1;
     }

#if 0
     int j;
     for (j = 0; j < realn_n_cigar; j++) {
          curr_oplen = realn_cigar[j] >> 4;
          curr_op = realn_cigar[j]&0xf;
          fprintf(stderr, "op:%d oplen:%d\n", curr_op, curr_oplen);
     }
#endif

     /* check if read was shifted */
     if (shift-(c->pos-lower) != 0) {
          LOG_VERBOSE("Read %s with shift of %d at original pos %s:%d\n", 
                      bam1_qname(b), shift-(c->pos-lower),
                      tmp->in->header->target_name[c->tid], c->pos);
          c->pos = c->pos + (shift - (c->pos - lower));
     }
     
	 if (reclip){
		 // check if first op or last op is I and replace with S
		 int curr_oplen_reclip = realn_cigar[0] >> 4;
		 int curr_op_reclip = realn_cigar[0]&0xf;

		 if (curr_op_reclip == BAM_CINS){
			curr_op_reclip = BAM_CSOFT_CLIP;
			realn_cigar[0] = curr_oplen_reclip <<4 | curr_op_reclip;
		}
		 curr_oplen_reclip = realn_cigar[realn_n_cigar-1] >> 4;
		 curr_op_reclip = realn_cigar[realn_n_cigar-1]&0xf;

		 if (curr_op_reclip == BAM_CINS){
			curr_op_reclip = BAM_CSOFT_CLIP;
			realn_cigar[realn_n_cigar-1] = curr_oplen_reclip <<4 | curr_op_reclip;
		}
	}
     replace_cigar(b, realn_n_cigar, realn_cigar);
     bam_write1(tmp->out, b);
     free(aln);
     free(realn_cigar);
     return 0;
}

static void usage()
{
     fprintf(stderr, "Usage: lofreq viterbi [options] in.bam\n");
     fprintf(stderr, "Options:\n");
     fprintf(stderr, "     -f | --ref FILE     Indexed reference fasta file [null]\n");
     fprintf(stderr, "     -k | --keepflags    Don't delete flags MC, MD, NM and A, which are all prone to change during realignment.\n");
     fprintf(stderr, "     -q | --defqual INT  Assume INT as quality for all bases with BQ2. Default (=-1) is to use median quality of bases in read.\n");
/* experimental. keep enabled but don't tell user about it */
#if 0
     fprintf(stderr, "     -r | --reclip       Reclip insertions and/or deletions on the beginning and end of read to soft clip\n");
     fprintf(stderr, "                         FILE HAS TO BE PREVIOUSLY UNCLIPPED!!!\n");
#endif
     fprintf(stderr, "     -o | --out FILE     Output BAM file [- = stdout = default]\n");
     fprintf(stderr, "          --verbose      Be verbose\n");
     fprintf(stderr, "\n");
     fprintf(stderr, "NOTE: Output BAM file will (likely) be unsorted (use samtools sort, e.g. lofreq viterbi ... | samtools sort -')\n");
}


int main_viterbi(int argc, char *argv[])
{
     tmpstruct_t tmp = {0};
     static int del_flag = 1;
     static int q2default = -1;
	 static int reclip = 0;
     char *bam_out = NULL;
     bam1_t *b = NULL;
 
     if (argc == 2) {
          usage();
          return 1;
     }

     while (1) {
          int c;

          static struct option long_options[] = {
               {"ref", required_argument, NULL, 'f'},
               {"verbose", no_argument, &verbose, 1},
               {"keepflags", no_argument, NULL, 'k'},
			   {"reclip",	no_argument, NULL, 'r'},
               {"out", required_argument, NULL, 'o'},
               {"defqual", required_argument, NULL, 'q'},
               {0,0,0,0}
          };
          
          static const char *long_opts_str = "rkf:q:o:";
          int long_option_index = 0;

          c = getopt_long(argc-1, argv+1, long_opts_str, long_options, &long_option_index);

          if (c == -1) {
               break;
          }
          switch (c){
          case 'f':
               if (! file_exists(optarg)) {
                    LOG_FATAL("Reference fasta file %s does not exist. Exiting...\n", optarg);
                    return 1;
               }
               tmp.fai = fai_load(optarg);	
               break;
          case 'k':
               del_flag = 0;
               break;
          case 'q':
               q2default = atoi(optarg);
               break;
		  case 'r':
				reclip = 1;
				break;
          case 'o':
               if (0 != strcmp(optarg, "-")) {
                    if (file_exists(optarg)) {
                         LOG_FATAL("Cowardly refusing to overwrite file '%s'. Exiting...\n", optarg);
                         return 1;
                    }
               }
               bam_out = strdup(optarg);
               break;
          case '?':
               LOG_FATAL("%s\n", "Unrecognized arguments found. Exiting\n");
               usage();
               break;
          default:
               break;
          }
     }


     if (! tmp.fai) {
          LOG_FATAL("%s\n", "Couldn't load reference fasta file\n");
          usage();
          return 1;
     }

     /* get bam file */
     if (1 != argc-optind-1){
          LOG_FATAL("%s\n", "Need exactly one BAM file as last argument\n");
          return 1;
     }
     if ((tmp.in = samopen((argv+optind+1)[0], "rb",0)) == 0){
          LOG_FATAL("Failed to open BAM file %s. Exiting...\n", (argv+optind+1)[0]);
          return 1;
     }

     if (!bam_out || bam_out[0] == '-') {
          tmp.out = bam_dopen(fileno(stdout), "w");
     } else {
          tmp.out = bam_open(bam_out, "w");
     }
     bam_header_write(tmp.out, tmp.in->header);
     
     b = bam_init1();
     tmp.tid = -1;
     tmp.ref = 0;
     while (samread(tmp.in, b) >= 0){
          fetch_func(b, &tmp, del_flag, q2default, reclip);
     }
     bam_destroy1(b);
     
     samclose(tmp.in);
     bam_close(tmp.out);
     if (tmp.ref)
          free(tmp.ref);
     fai_destroy(tmp.fai);
     free(bam_out);

     LOG_VERBOSE("%s\n", "NOTE: Output BAM file will be unsorted (use samtools sort, e.g. samtools sort -')");

     return 0;
}
