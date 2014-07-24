#ifndef VITERBI_H
#define VITERBI_H
int left_align_indels(char *sref, char *squery, int slen, char *res);
int viterbi(char *ref, char *query, char *bqual, char *aln, int quality);
int viterbi_test();
#endif
