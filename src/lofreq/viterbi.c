#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "viterbi.h"

#define SANGERQUAL_TO_PROB(qual) (pow(10.0, -0.1*((int)qual - 33)))

int argmax(double *array, int alen) {
     int i;
     int curr_i = 0;
     double curr_max = INT_MIN;
     for (i = 0; i < alen; i++) {
          if (array[i] > curr_max) {
               curr_i = i;
               curr_max = array[i];
          }
     }
     return curr_i;
}

int left_align_indels(char *sref, char *squery, int slen, char *new_state_seq) {

     char ref[slen];
     char query[slen];
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

int viterbi(char *ref, char *query, char *bqual, char *aln)
{
     int qlen = strlen(query)+1;
     int rlen = strlen(ref)+1;
          
     // Define transition probabilities
     // FIXME: define globally to speed up
     double alpha = 0.00001;
     double beta = 0.4;
     double L = (double)rlen;
     double gamma = 1/(2.*L);

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
     
     double ep_ins = log10(.25); // Insertion emission probability

     // Initialize

     double V_start[qlen];
     double V_match[rlen][qlen];
     double V_ins[rlen][qlen];
     double V_del[rlen][qlen];

     memset(V_match, 0, sizeof(V_match));
     memset(V_ins, 0, sizeof(V_ins));
     memset(V_del, 0, sizeof(V_del));
     
     int i, k;
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

     char ptr_match[rlen][qlen];
     char ptr_ins[rlen][qlen];
     char ptr_del[rlen][qlen];

     memset(ptr_match, 0, sizeof(ptr_match));
     memset(ptr_ins, 0, sizeof(ptr_ins));
     memset(ptr_del, 0, sizeof(ptr_del));

     // Recursion

     double bp;
     for (i = 1; i < qlen; i++) {
          
          // Define emission probabilities
           
          bp = SANGERQUAL_TO_PROB(bqual[i-1]);
          
          double ep_match = log10(1-bp);
          double ep_match_not = log10(bp/3.);

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
               index = argmax(mterms, 4);
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
               index = argmax(iterms, 3);
               ptr_ins[k][i] = "SMI"[index];
               V_ins[k][i] = ep_ins + iterms[index];

               // V_Dk(i) = max( M_k-1(i) + log(a_(M_k-1,D_k)),
               //                D_k-1(i) + log(a_(D_k-1,D_k)) )
               double dterms[2] = {V_match[k-1][i] + tp[0][2],
                                  V_del[k-1][i] + tp[2][2]};
               index = argmax(dterms, 2);
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
     return k;
}

int viterbi_test()
{

     fprintf(stderr, "Testing viterbi realignment...\n");
     viterbi("CCATATGG", "CCATGG", "??????", NULL);

     fprintf(stderr, "Testing left-alignment of indels...\n");
     left_align_indels("CCATATGG", "CCAT**GG", 8, 0);
     left_align_indels("CCAT**GG", "CCATATGG", 8, 0);
     left_align_indels("CCATATGG*CC", "CCAT**GGGCC", 11, 0);

     return 0;
}
