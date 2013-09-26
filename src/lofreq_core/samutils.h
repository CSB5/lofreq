#ifndef SAMUTILS_H
#define SAMUTILS_H


typedef enum {
        OP_MATCH,
        OP_MISMATCH,
        OP_INS,
        OP_DEL,
        NUM_OP_CATS,
} op_cat_t;

#define STR(name) # name

static char *op_cat_str[] = {
    STR(OP_MATCH),
    STR(OP_MISMATCH),
    STR(OP_INS),
    STR(OP_DEL),
    STR(NUM_OP_CATS)
};


char *
cigar_str_from_bam(const bam1_t *b);

int
count_matches(int *counts, const bam1_t *b, const char *ref, int min_bq);


#endif
