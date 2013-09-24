#ifndef SAMUTILS_H
#define SAMUTILS_H

char *
cigar_from_bam(const bam1_t *b);

#define MATCH_COUNT_IDX 0
#define MISMATCH_COUNT_IDX 1
#define INS_COUNT_IDX 2
#define DEL_COUNT_IDX 3
#define MAX_COUNT_IDX 4

int
count_matches(int *counts, const bam1_t *b, const char *ref, int min_bq);


#endif
