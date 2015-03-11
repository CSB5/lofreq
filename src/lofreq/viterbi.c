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

#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "viterbi.h"
#include "utils.h"

#define PHRED_TO_SANGERQUAL(i) ((char)(i)+33)
#define SANGERQUAL_TO_PHRED(c) ((int)(c)-33)
#define SANGERQUAL_TO_PROB(c)  (pow(10.0, -0.1*SANGERQUAL_TO_PHRED(c)))

const int BQ2_DEFAULT = 20;

#ifdef PACBIO_REALN
static int pacbio_msg_printed = 0;
#endif

int left_align_indels(char *sref, char *squery, int slen, char *new_state_seq) {

     char ref[slen+1];
     char query[slen+1];
     strcpy(ref, sref);
     strcpy(query, squery);
     
     int i = 0;
     // FIXME: can be further optimized
     while (i < slen-1) {
          if (ref[i] != '*' && query[i] != '*') {
               if (ref[i+1] == '*') {
                    int ilen = 0;
                    while (ref[i+1+ilen] == '*') { ilen++; }
                    if (query[i+ilen] == ref[i]) {
                         ref[i+ilen] = ref[i];
                         ref[i] = '*';
                         i--;
                         continue;
                    }
               } else if (query[i+1] == '*') {
                    int dlen = 0;
                    while (query[i+1+dlen] == '*') { dlen++; }
                    if (query[i] == ref[i+dlen]) {
                         query[i+dlen] = query[i];
                         query[i] = '*';
                         i--;
                         continue;
                    }
               }
          }
          i++;
     }

     char state_seq[slen+1];
     for (i = 0; i < slen; i++) {
          if (ref[i] == '*') { state_seq[i] = 'I'; }
          else if (query[i] == '*') { state_seq[i] = 'D'; }
          else { state_seq[i] = 'M'; }
     }
     state_seq[i] = '\0';
     //fprintf(stderr, "ref:%s, query:%s, state_seq:%s\n", ref, query, state_seq);
     
     if (new_state_seq) {
          strcpy(new_state_seq, state_seq);
     }
     
     return 0;
}

/* bqual is the base quality phred score representation as string. so use SANGERQUAL_TO_PROB for conversion */
int viterbi(char *ref, char *query, char *bqual, char *aln, int quality)
{
     //printf("inside viterbi\n");
     int qlen = strlen(query)+1;
     int rlen = strlen(ref)+1;
     
     double *V_start;
     double **V_match;
     double **V_ins;
     double **V_del;
          
     char **ptr_match;
     char **ptr_ins;
     char **ptr_del;

     // Define transition probabilities
     // FIXME: define globally to speed up
#ifdef PACBIO_REALN
     double alpha = 0.1;
     if (! pacbio_msg_printed) {
          fprintf(stderr, "WARN(%s|%s): Using pacbio viterbi params\n", __FILE__, __FUNCTION__);
          pacbio_msg_printed = 1;
     }
#else
     double alpha = 0.00001;
#endif
     double beta = 0.4;

     double L = (double)rlen;
     double gamma = 1/(2.*L);
     int i, k;
     double ep_ins = log10(.25); // Insertion emission probability
     double tp[5][5] = {{0}};

     tp[0][0] = log10((1 - 2*alpha)*(1 - gamma)); // M->M
     tp[0][1] = log10(alpha*(1 - gamma)); // M->I
     tp[0][2] = log10(alpha*(1 - gamma)); // M->D
     tp[0][4] = log10(gamma); // M->E
     tp[1][0] = log10((1 - beta)*(1 - gamma)); // I->M
     tp[1][1] = log10(beta*(1 - gamma)); // I->I
     tp[1][4] = log10(gamma); // I->E
     tp[2][0] = log10(1- beta); // D->M
     tp[2][2] = log10(beta); // D->D
     tp[3][0] = log10((1 - alpha)/L); // S->M
     tp[3][1] = log10(alpha/L); // S->I
     

     // Initialize
     V_start = malloc(qlen * sizeof(double));
     V_match = malloc(rlen * sizeof(double*));
     V_ins = malloc(rlen * sizeof(double*));
     V_del = malloc(rlen * sizeof(double*));
     for (i = 0; i < rlen; i++) { 
          V_match[i] = calloc(qlen, sizeof(double));
          V_ins[i] = calloc(qlen, sizeof(double));
          V_del[i] = calloc(qlen, sizeof(double));
     }

     for (i = 0; i < qlen; i++) { 
          V_start[i] = INT_MIN; 
     }
     for (k = 0; k < rlen; k++) { 
          V_match[k][0] = INT_MIN; 
          V_ins[k][0] = INT_MIN;
          V_del[k][0] = INT_MIN;
     }
     for (i = 0; i < qlen; i++) { 
          V_match[0][i] = INT_MIN; 
          V_ins[0][i] = INT_MIN;
          V_del[0][i] = INT_MIN;
     }
     V_start[0] = 0;

     ptr_match = malloc(rlen * sizeof(char*));
     ptr_ins = malloc(rlen * sizeof(char*));
     ptr_del = malloc(rlen * sizeof(char*));
     for (i=0; i<rlen; i++) {
          ptr_match[i] = calloc(qlen, sizeof(char));
          ptr_ins[i] = calloc(qlen, sizeof(char));
          ptr_del[i] = calloc(qlen, sizeof(char));
     }    
    
     // Recursion
     double bp;
     for (i = 1; i < qlen; i++) {
          double ep_match;
          double ep_match_not;

          // Define emission probabilities
		  if ( SANGERQUAL_TO_PHRED(bqual[i-1]) == 2) {
               bp = SANGERQUAL_TO_PROB(PHRED_TO_SANGERQUAL(quality));
		  } else {
               bp = SANGERQUAL_TO_PROB(bqual[i-1]);
		  }
          ep_match = log10(1-bp);
          ep_match_not = log10(bp/3.);
                    
          for (k = 1; k < rlen; k++) {
               int index;
               
               // V_Mk(i) = log(e_Mk(x_i)) + max( S_0(i-1) + log(a_(S_0,M_k)),
               //                                 M_k-1(i-1) + log(a_(M_k-1,M_k)),
               //                                 I_k-1(i-1) + log(a_(I_k-1,M_k)),
               //                                 D_k-1(i-1) + log(a_(D_k-1,M_k)) )
               double mterms[4] = {V_start[i-1] + tp[3][0],
                                   V_match[k-1][i-1] + tp[0][0],
                                   V_ins[k-1][i-1] + tp[1][0],
                                   V_del[k-1][i-1] + tp[2][0]};
               index = argmax_d(mterms, 4);
               ptr_match[k][i] = "SMID"[index];
               if (query[i-1] == ref[k-1]) {
                    V_match[k][i] = ep_match + mterms[index];
               } else {
                    V_match[k][i] = ep_match_not + mterms[index];
               }
               
               // V_Ik(i) = log(e_Ik(x_i)) + max( S_0(i-1) + log(a_(S_0,I_k)),
               //                                 M_k(i-1) + log(a_(M_k,I_k)),
               //                                 I_k(i-1) + log(a_(I_k,I_k)) )
               
               double iterms[3] = {V_start[i-1] + tp[3][1],
                                   V_match[k][i-1] + tp[0][1],
                                   V_ins[k][i-1] + tp[1][1]};
               index = argmax_d(iterms, 3);
               ptr_ins[k][i] = "SMI"[index];
               V_ins[k][i] = ep_ins + iterms[index];
               
               // V_Dk(i) = max( M_k-1(i) + log(a_(M_k-1,D_k)),
               //                D_k-1(i) + log(a_(D_k-1,D_k)) )
               double dterms[2] = {V_match[k-1][i] + tp[0][2],
                                   V_del[k-1][i] + tp[2][2]};
               index = argmax_d(dterms, 2);
               ptr_del[k][i] = "MD"[index];
               V_del[k][i] = dterms[index];
               
               //fprintf(stderr, "k:%d, i:%d, %f, %f, %f\n", k, i, 
               //                 V_match[k][i], V_ins[k][i], V_del[k][i]);
          }
     }
     
     // Termination
     // max[M_L(N), I_L(N), D_L(N)]
     char end_state = '!';
     double best_score = INT_MIN;
     int best_index = 0;
     for (k = 0; k < rlen; k++) {
          if (V_match[k][qlen-1] > best_score) {
               end_state = 'M';
               best_score = V_match[k][qlen-1];
               best_index = k;
          }
          if (V_ins[k][qlen-1] > best_score) {
               end_state = 'I';
               best_score = V_ins[k][qlen-1];
               best_index = k;
          }
     }
     //fprintf(stderr, "ended on %c, best_score is %f, best_index is %d\n",
     //     end_state, best_score, best_index);
     for (i = 0; i < rlen; i++) { 
          free(V_match[i]);
          free(V_ins[i]);
          free(V_del[i]);
     }
     free(V_match);
     free(V_ins);
     free(V_del);
     free(V_start);

     // Trace-back
     i = qlen - 1;
     k = best_index;
     int maxslen = qlen+rlen;
     char current_ptr = end_state;
     char tmp_state_seq[maxslen], tmp_ref[maxslen], tmp_query[maxslen];
     tmp_state_seq[qlen+rlen-1] = tmp_ref[qlen+rlen-1] = tmp_query[qlen+rlen-1] = '\0';
     int si = qlen+rlen-2;
     
     while (i != 0 && k != 0) {
          tmp_state_seq[si] = current_ptr;
          if (current_ptr == 'S') {
               break;
          } else if (current_ptr == 'M') {
               tmp_ref[si] = ref[k-1];
               tmp_query[si] = query[i-1];
               current_ptr = ptr_match[k][i];
               i -= 1;
               k -= 1;
          } else if (current_ptr == 'I') {
               tmp_ref[si] = '*';
               tmp_query[si] = query[i-1];
               current_ptr = ptr_ins[k][i];
               i -= 1;
          } else if (current_ptr == 'D') {
               tmp_ref[si] = ref[k-1];
               tmp_query[si] = '*';
               current_ptr = ptr_del[k][i];
               k -= 1;
          } else {
               return -1;
          }
          si--;
     }
     for (i=0; i<rlen; i++) {
          free(ptr_match[i]);
          free(ptr_ins[i]);
          free(ptr_del[i]);
     }
     free(ptr_match);
     free(ptr_ins);
     free(ptr_del);


     {
          char *state_seq = tmp_state_seq+si+1;
          char *new_ref = tmp_ref+si+1;
          char *new_query = tmp_query+si+1;
          //fprintf(stderr, "ref:%s, query:%s, state_seq:%s\n", ref+1, query+1, state_seq);
          int state_seq_len = strlen(state_seq);
          char *new_state_seq = malloc(sizeof(char)*(state_seq_len+1));
          left_align_indels(new_ref, new_query, state_seq_len, new_state_seq);
          
          if (aln) {
               strcpy(aln, new_state_seq);
          }
     
          free(new_state_seq);
     }

    return k;

}

int viterbi_test()
{

     fprintf(stderr, "Testing viterbi realignment...\n");
     viterbi("CCATATGG", "CCATGG", "??????", NULL, 20);

     fprintf(stderr, "Testing left-alignment of indels...\n");
     left_align_indels("CCATATGG", "CCAT**GG", 8, 0);
     left_align_indels("CCAT**GG", "CCATATGG", 8, 0);
     left_align_indels("CCATATGG*CC", "CCAT**GGGCC", 11, 0);

     return 0;
}

//test to figure out quality calculation
#ifdef VITERBI_MAIN
int main(int argc, char *argv[])
{
  int i;
  char *bq = "#5?";
  for (i=0; i<strlen(bq); i++) {
    char c = bq[i];
    printf("#%d: %c=%g=%d (cast %d)\n", i+1, c, SANGERQUAL_TO_PROB(c), SANGERQUAL_TO_PHRED(c), (int)c);
  }
  return 0;
}
#endif
