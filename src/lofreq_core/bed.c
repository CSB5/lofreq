/* -*- c-file-style: "k&r"; indent-tabs-mode: nil; -*- */

/* FIXME can we replace all of this libbam's bed_read() etc.
 * FIXME doesn't behave as expected if regions overlap
 */

#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <assert.h>

#include "log.h"
#include "bed.h"

#define BUF_SIZE 1024


void
free_bed(bed_t *bed)
{
    int i;
    for (i=0; i<bed->nregions; i++) {
        free(bed->region[i].chrom);
    }
    free(bed->region);
}


/* returns -1 on error. bed_pos_t members allocated here. have to be
 * freed by caller.
 */
int
parse_bed(bed_t *bed, const char *bed_file) {
    char line[BUF_SIZE];
    FILE *fh;

    assert(0==bed->nregions);
    bed->nregions = 0;
    bed->region = NULL;

    fh = fopen(bed_file, "r");
    if (NULL == fh) {
        LOG_ERROR("Couldn't open bed-file %s\n", bed_file);
        return -1;
    }

    while (NULL != fgets(line, BUF_SIZE, fh)) {
        char chrom[BUF_SIZE];
        long int start, end;
        if (line[0]=='#') {
            continue;
        }
        if (1 == strlen(line)) {
            LOG_WARN("Skippping empty line in bed-file %s", bed_file);
            continue;
        }

        /* this works with any number of tabs and white-spaces */
        if (3 != sscanf(line, "%s %ld %ld \n", chrom, &start, &end)) {
            LOG_FATAL("Couldn't parse the following line"
                      " from bed-file %s: %s", bed_file, line);
            fclose(fh);
            return -1;
        }
        if (! (end>start)) {
            LOG_WARN("ignoring the following bogus line"
                      " in bed-file %s: %s", bed_file, line);
            continue;
        }

        bed->nregions += 1;
        bed->region = realloc(bed->region, bed->nregions*sizeof(region_t));
        bed->region[bed->nregions-1].chrom = strdup(chrom);
        bed->region[bed->nregions-1].start = start;
        bed->region[bed->nregions-1].end = end;
    }
    fclose(fh);
    return bed->nregions;
}
/* parse_bed */


void
dump_bed(bed_t *bed)
{
    int i;
    for (i=0; i<bed->nregions; i++) {
        fprintf(stderr, "%d: %s %d %d\n", i,
                bed->region[i].chrom,
                bed->region[i].start,
                bed->region[i].end);
    }
}


/* returns -1 on error 
 */
long long int
bed_pos_sum(bed_t *bed) {
    long long int sum = 0;
    int i;
    for (i=0; i<bed->nregions; i++) {
        if (sum > LLONG_MAX - (bed->region[i].end - bed->region[i].start)) {
            LOG_FATAL("%s\n", "count overflow!");
            return -1;
        }
        sum += (bed->region[i].end - bed->region[i].start);
    }
    return sum;
}
/* bed_pos_sum */

