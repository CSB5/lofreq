/* -*- c-file-style: "k&r"; indent-tabs-mode: nil; -*- */
/*********************************************************************
 *
 * FIXME update license
 *
 *********************************************************************/
#include <ctype.h>
#include <stdio.h>
#include <getopt.h>

#include "faidx.h"
#include "sam.h" 
#include "log.h"
#include "utils.h"
#include "lofreq_indel_quality.h"


char DINDELQ[] = "!MMMLKEC@=<;:988776"; /* 1-based 18 */
char DINDELQ2[] = "!CCCBA;963210/----,"; /* *10 */


typedef struct {
     samfile_t *in;
     bamFile out;
     uint8_t defq[1000];
} tmpstruct_t_default;


typedef struct {
     samfile_t *in;
     bamFile out;
     faidx_t *fai;
     int *hpcount;
     int rlen;
     uint32_t tid;
} tmpstruct_t_dindel;


#define ENCODE_Q(q) (uint8_t)(q < 33 ? '!' : (q > 126 ? '~' : q))


static int default_fetch_func(bam1_t *b, void *data)
{
     tmpstruct_t_default *tmp = (tmpstruct_t_default*)data;
     bam1_core_t *c = &b->core;

     uint8_t indelq[c->l_qseq+1];

     memcpy(&indelq, &tmp->defq, c->l_qseq);
     indelq[c->l_qseq] = '\0';

     uint8_t *to_delete;
     to_delete = bam_aux_get(b, "BI");
     if (to_delete) bam_aux_del(b, to_delete);
     bam_aux_append(b, "BI", 'Z', c->l_qseq+1, indelq);


     to_delete = bam_aux_get(b, "BD");
     if (to_delete) bam_aux_del(b, to_delete);
     bam_aux_append(b, "BD", 'Z', c->l_qseq+1, indelq);

     bam_write1(tmp->out, b);
     return 0;
}


/* Stores an array of ints that corresponds to the length of the
 * homopolymer at the start of each homopolymer*/
int find_homopolymers(char *query, int *count, int qlen)
{
     int i, j;
     int curr_i = 0;
     int curr_count = 1;
     for (i = 1; i < qlen; i++) {
          if (query[i] == query[curr_i]) {
                  curr_count += 1; /* keep incrementing count if in homopolymer region */
          } else {
                  count[curr_i] = curr_count; /* record length of homopolymer region */
               for (j = curr_i+1; j < i; j++) {
                       count[j] = 1; /* all other positions get a count of 1 */
               }
               curr_i = i;
               curr_count = 1;
          }
     }
     if (curr_i < i) { /* take care of edge case at the end of the read */
          count[curr_i] = curr_count;
          for (j = curr_i+1; j < i; j++) {
               count[j] = 1;
          }
     }
     return 0;
}


static int dindel_fetch_func(bam1_t *b, void *data)
{
     tmpstruct_t_dindel *tmp = (tmpstruct_t_dindel*)data;
     bam1_core_t *c = &b->core;
     int rlen;

     /* skip all reads that are not properly paired */
     if (! (c->flag & BAM_FPROPER_PAIR)) {
             /* fprintf(stderr, "skipping read: %s at pos %d\n", bam1_qname(b), c->pos); */
          return 0;
     }

     /* get the reference sequence and compute homopolymer array */
     if (tmp->tid != c->tid) {
             /*fprintf(stderr, "fetching reference sequence %s\n",
               tmp->in->header->target_name[c->tid]); */
          char *ref = fai_fetch(tmp->fai, tmp->in->header->target_name[c->tid], &rlen);
          int rlen = strlen(ref);
          tmp->tid = c->tid;
          if (tmp->hpcount) free(tmp->hpcount);
          tmp->hpcount = (int*)malloc(rlen*sizeof(int));
          find_homopolymers(ref, tmp->hpcount, rlen);
          free(ref);
          tmp->rlen = rlen;
          /* fprintf(stderr, "fetched reference sequence\n");*/
     }

     /* parse the cigar string */
     uint32_t *cigar = bam1_cigar(b);
     uint8_t indelq[c->l_qseq+1];
     /* fprintf(stderr, "l_qseq:%d\n", c->l_qseq); */
     int i;
     int x = c->pos; /* coordinate on reference */
     int y = 0; /* coordinate on query */
     for (i = 0; i < c->n_cigar; ++i) {
          int j, oplen = cigar[i]>>4, op = cigar[i]&0xf;
          if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF) {
               for (j = 0; j < oplen; j++) {
                       /*fprintf(stderr, "query:%d, ref:%d, count:%d\n", 
                         y, x, tmp->hpcount[x+1]); */
                    indelq[y] = (x > tmp->rlen-2) ? DINDELQ[0] : (tmp->hpcount[x+1]>18? 
                         DINDELQ[0] : DINDELQ[tmp->hpcount[x+1]]);
                    x++; 
                    y++;
               }
          } else if (op == BAM_CHARD_CLIP) { /* do nothing */
          } else if (op == BAM_CDEL) {
               x += oplen;
          } else if (op == BAM_CINS || op == BAM_CSOFT_CLIP) { 
               for (j = 0; j < oplen; j++) {
                       /* fprintf(stderr, "query:%d, ref:%d\n", y, x); */
                    indelq[y] = DINDELQ[0];
                    y++;
               }
          } else {
               exit(1);
          }
     }
     indelq[y] = '\0';

     bam_aux_append(b, "BI", 'Z', c->l_qseq+1, indelq);
     bam_aux_append(b, "BD", 'Z', c->l_qseq+1, indelq);

     bam_write1(tmp->out, b);
     return 0;
}


int add_uniform(const char *bam_in, const char *bam_out, const int qual)
{
	tmpstruct_t_default tmp;
    uint8_t q = ENCODE_Q(qual+33);
    int i;
    bam1_t *b = NULL;
    int count = 0;

	if ((tmp.in = samopen(bam_in, "rb", 0)) == 0) {
         LOG_FATAL("Failed to open BAM file %s\n", bam_in);
         return 1;
    }

    for (i = 0; i < 999; i++) {
         tmp.defq[i] = q;
    }
    tmp.defq[i] = '\0';

    if (!bam_out || bam_out[0] == '-') {
         tmp.out = bam_dopen(fileno(stdout), "w");
    } else {
         tmp.out = bam_open(bam_out, "w");
    }
    bam_header_write(tmp.out, tmp.in->header);
    
    b = bam_init1();
    while (samread(tmp.in, b) >= 0) {
         count++;
         default_fetch_func(b, &tmp); 
    }
    bam_destroy1(b);
    
    samclose(tmp.in);
    bam_close(tmp.out);
    LOG_VERBOSE("Processed %d reads\n", count);
    return 0;
}


int add_dindel(const char *bam_in, const char *bam_out, const char *ref)
{
	tmpstruct_t_dindel tmp;
    int count = 0;
    bam1_t *b = NULL;

	if ((tmp.in = samopen(bam_in, "rb", 0)) == 0) {
         LOG_FATAL("Failed to open BAM file %s\n", bam_in);
             return 1;
        }
    if ((tmp.fai = fai_load(ref)) == 0) {
         LOG_FATAL("Failed to open reference file %s\n", ref);
         return 1;
    }
    /*warn_old_fai(ref);*/

    if (!bam_out || bam_out[0] == '-') {
         tmp.out = bam_dopen(fileno(stdout), "w");
    } else {
         tmp.out = bam_open(bam_out, "w");
    }
    bam_header_write(tmp.out, tmp.in->header);
    
    b = bam_init1();
    tmp.tid = -1;
    tmp.hpcount = 0;
    tmp.rlen = 0;
    while (samread(tmp.in, b) >= 0) {
         count++;
         dindel_fetch_func(b, &tmp); 
    }
    bam_destroy1(b);
    
    if (tmp.hpcount) free(tmp.hpcount);
    samclose(tmp.in);
    bam_close(tmp.out);
    fai_destroy(tmp.fai);
	LOG_VERBOSE("Processed %d reads\n", count);
	return 0;
}


static void
usage()
{
     const char *myname = "lofreq indel_quality";
     fprintf(stderr,
             "%s: Insert indel qualities into BAM file (required for indel predictions)." \
             "\n\n"  \
             "The preferred way of doing this is via GATK's BQSR!" \
             " If that's not possible, use this subcommand.\n"  \
             "The command has two modes: 'uniform' and 'dindel':\n" \
             "- 'uniform' will assign a given value uniformly, whereas\n"  \
             "- 'dindel' will insert indel qualities based on Dindel (PMID 20980555).\n\n", myname);

     fprintf(stderr,"Usage: %s [options] in.bam\n\n", myname);
     fprintf(stderr,"Options:\n");
     fprintf(stderr, "  -u | --uniform INT  Add this indel quality uniformly to all bases");
     fprintf(stderr, " (clashes with --dindel)\n");
     fprintf(stderr, "       --dindel       Add Dindel's indel qualities");
     fprintf(stderr, " (clashes with -u; needs --ref)\n");
     fprintf(stderr, "  -f | --ref          Reference sequence used for mapping");
     fprintf(stderr, "  -o | --out FILE     Output BAM file [- = stdout = default]\n");
     fprintf(stderr, "\n\n");
     fprintf(stderr, "Do not realign your BAM file after this!\n");
     fprintf(stderr, "\n");
}



int main_indel_quality(int argc, char *argv[])
{
     char *bam_in = NULL;
     char *bam_out = NULL; /* - == stdout */
     char *ref = NULL;
     int c;
     static int dindel = 0;
     int uniform = -1;
     while (1) {
          static struct option long_opts[] = {
               /* see usage sync */
               {"help", no_argument, NULL, 'h'},
               {"verbose", no_argument, &verbose, 1},
               {"debug", no_argument, &debug, 1},
               {"dindel", no_argument, &dindel, 1},
               {"out", required_argument, NULL, 'o'},
               {"uniform", required_argument, NULL, 'u'},
               {"ref", required_argument, NULL, 'f'},
               {0, 0, 0, 0} /* sentinel */
          };
          
          /* keep in sync with long_opts and usage */
          static const char *long_opts_str = "hu:f:o:";
     
          /* getopt_long stores the option index here. */
          int long_opts_index = 0;
          c = getopt_long(argc-1, argv+1, /* skipping 'lofreq', just leaving 'command', i.e. call */
                          long_opts_str, long_opts, & long_opts_index);
          if (c == -1) {
               break;
          }
          switch (c) {
               /* keep in sync with long_opts etc */
          case 'h':
               usage();
               return 0;
          case 'u':
               uniform = atoi(optarg);
               break;
          case 'f':
               if (! file_exists(optarg)) {
                    LOG_FATAL("Reference fasta file '%s' does not exist. Exiting...\n", optarg);
                    return 1;
               }
              ref = strdup(optarg);
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
               LOG_FATAL("%s\n", "unrecognized arguments found. Exiting...\n");
               return 1;
          default:
               break;
          }
     }
     if (1 != argc - optind - 1) {
          fprintf(stderr, "Need exactly one BAM file as last argument\n");
          return 1;
     }
     bam_in = (argv + optind + 1)[0];
     if (! file_exists(bam_in)) {
          LOG_FATAL("BAM file %s does not exist. Exiting...\n", bam_in);
          return -1;
     }

     if (! bam_out) {
          bam_out = malloc(2 * sizeof(char));
          strcpy(bam_out, "-");
     }

     if (uniform != -1) {
          if (dindel) {
               LOG_FATAL("%s\n", "Can't insert both, uniform and dindel qualities");
               return -1;
          }
          return add_uniform(bam_in, bam_out, uniform);          

     } else if (dindel) {
          if (! ref) {
               LOG_FATAL("%s\n", "Need reference for Dindel model");
               return -1;
          }
          return add_dindel(bam_in, bam_out, ref);          

     } else {
          LOG_FATAL("%s\n", "Please specify either dindel or uniform mode");
          return -1;
     }
     free(ref);
     free(bam_out);
     return 0;
}
