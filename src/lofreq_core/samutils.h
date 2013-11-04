#ifndef SAMUTILS_H
#define SAMUTILS_H

/* FIXME should become shared const MAX_READ_LEN. values <10k might not be enough for pacbio */
#define MAX_READ_LEN 8192


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
count_cigar_ops(int *counts, int **quals,
                const bam1_t *b, const char *ref, int min_bq,
                char *target);


#ifdef USE_MAPERRPROF

typedef struct {
     int num_targets; /* bam_header->n_targets */
     int *prop_len; /* one prop length per target: index is tid */
     double **props; /* one prop array per target: index is tid */
} alnerrprof_t;


void
normalize_alnerrprof(alnerrprof_t *alnerrprof);

int
parse_alnerrprof_statsfile(alnerrprof_t *alnerrprof, const char *path, bam_header_t *bam_header);

void
calc_read_alnerrprof(double *alnerrprof, unsigned long int *used_pos, 
                        const bam1_t *b, const char *ref);

void
write_alnerrprof_stats(char *target_name, unsigned long int *alnerrprof_usedpos, 
                    double *alnerrprof, int max_obs_read_len, FILE *out);

void
free_alnerrprof(alnerrprof_t *alnerrprof);

#endif


#endif
