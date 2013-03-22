/* -*- c-file-style: "k&r" -*- */
#ifndef BED_H
#define BED_H


typedef struct {
     char *chrom;
     int start; /* zero-offset */
     int end; /* excluding */
} region_t;

typedef struct {
     region_t *region;
     int nregions;
} bed_t;

int parse_bed(bed_t *b, const char *f);

long long int bed_pos_sum(bed_t *b);

void dump_bed(bed_t *bed);

void free_bed(bed_t *bed);

#endif
