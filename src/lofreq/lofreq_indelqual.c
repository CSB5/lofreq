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
#include <string.h>
#include <getopt.h>

#include "htslib/faidx.h"
#include "sam.h" 
#include "log.h"
#include "utils.h"
#include "defaults.h"
#include "lofreq_indelqual.h"


char DINDELQ[] = "!MMMLKEC@=<;:988776"; /* 1-based 18 */
char DINDELQ2[] = "!CCCBA;963210/----,"; /* *10 */


typedef struct {
     samfile_t *in;
     bamFile out;
     int iq;
     int dq;
} data_t_uniform;


typedef struct {
     samfile_t *in;
     bamFile out;
     faidx_t *fai;
     int *hpcount;
     int rlen;
     uint32_t tid;
} data_t_dindel;


#define ENCODE_Q(q) (uint8_t)(q < 33 ? '!' : (q > 126 ? '~' : q))


static int uniform_fetch_func(bam1_t *b, void *data)
{
     uint8_t *to_delete;
     data_t_uniform *tmp = (data_t_uniform*)data;
     bam1_core_t *c = &b->core;
     char *iq;
     char *dq;

     iq = malloc((c->l_qseq+1) * sizeof(char));
     memset(iq, tmp->iq, c->l_qseq);
     iq[c->l_qseq] = '\0';

     to_delete = bam_aux_get(b, BI_TAG);
     if (to_delete) {
          bam_aux_del(b, to_delete);
     }
     bam_aux_append(b, BI_TAG, 'Z', c->l_qseq+1, (uint8_t*) iq);


     dq = malloc((c->l_qseq+1) * sizeof(char));
     memset(dq, tmp->dq, c->l_qseq);
     dq[c->l_qseq] = '\0';

     to_delete = bam_aux_get(b, BD_TAG);
     if (to_delete) {
          bam_aux_del(b, to_delete);
     }
     bam_aux_append(b, BD_TAG, 'Z', c->l_qseq+1, (uint8_t*) dq);

     bam_write1(tmp->out, b);

     free(iq);
     free(dq);

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
     data_t_dindel *tmp = (data_t_dindel*)data;
     bam1_core_t *c = &b->core;
     int rlen;
     uint8_t *to_delete;

     /* don't change reads failing default mask: BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP */
     if (c->flag & BAM_DEF_MASK) {
          /* fprintf(stderr, "skipping read: %s at pos %d\n", bam1_qname(b), c->pos); */
          bam_write1(tmp->out, b);
          return 0;
     }

     /* get the reference sequence and compute homopolymer array */
     if (tmp->tid != c->tid) {
             /*fprintf(stderr, "fetching reference sequence %s\n",
               tmp->in->header->target_name[c->tid]); */
          char *ref = fai_fetch(tmp->fai, tmp->in->header->target_name[c->tid], &rlen);
          strtoupper(ref);/* safeguard */
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
                    /* FIXME clang complains: The left operand of '>' is a garbage value */
                    indelq[y] = (x > tmp->rlen-2) ? DINDELQ[0] : (tmp->hpcount[x+1]>18 ?
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
               LOG_FATAL("unknown op %d for read %s\n", op, bam1_qname(b));/* FIXME skip? seen this somewhere else properly handled */
               exit(1);
          }
     }
     indelq[y] = '\0';

     to_delete = bam_aux_get(b, BI_TAG);
     if (to_delete) {
          bam_aux_del(b, to_delete);
     }
     bam_aux_append(b, BI_TAG, 'Z', c->l_qseq+1, indelq);

     to_delete = bam_aux_get(b, BD_TAG);
     if (to_delete) {
          bam_aux_del(b, to_delete);
     }
     bam_aux_append(b, BD_TAG, 'Z', c->l_qseq+1, indelq);

     bam_write1(tmp->out, b);
     return 0;
}


int add_uniform(const char *bam_in, const char *bam_out,
                const int ins_qual, const int del_qual)
{
	data_t_uniform tmp;
    uint8_t iq = ENCODE_Q(ins_qual+33);
    uint8_t dq = ENCODE_Q(del_qual+33);
    bam1_t *b = NULL;
    int count = 0;

	if ((tmp.in = samopen(bam_in, "rb", 0)) == 0) {
         LOG_FATAL("Failed to open BAM file %s\n", bam_in);
         return 1;
    }

    tmp.iq = iq;
    tmp.dq = dq;

    if (!bam_out || bam_out[0] == '-') {
         tmp.out = bam_dopen(fileno(stdout), "w");
    } else {
         tmp.out = bam_open(bam_out, "w");
    }
    bam_header_write(tmp.out, tmp.in->header);
    
    b = bam_init1();
    while (samread(tmp.in, b) >= 0) {
         count++;
         uniform_fetch_func(b, &tmp); 
    }
    bam_destroy1(b);
    
    samclose(tmp.in);
    bam_close(tmp.out);
    LOG_VERBOSE("Processed %d reads\n", count);
    return 0;
}


int add_dindel(const char *bam_in, const char *bam_out, const char *ref)
{
	data_t_dindel tmp;
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
     const char *myname = "lofreq indelqual";
     fprintf(stderr, "%s: Insert indel qualities into BAM file (required for indel predictions)\n\n", myname);
     fprintf(stderr, "Usage: %s [options] in.bam\n", myname);
     fprintf(stderr,"Options:\n");
     fprintf(stderr, "  -u | --uniform INT[,INT]  Add this indel quality uniformly to all bases.\n");
     fprintf(stderr, "                            Use two comma separated values to specify\n");
     fprintf(stderr, "                            insertion and deletion quality separately.\n");
     fprintf(stderr, "                            (clashes with --dindel)\n");
     fprintf(stderr, "       --dindel             Add Dindel's indel qualities (Illumina specific)\n");
     fprintf(stderr, "                            (clashes with -u; needs --ref)\n");
     fprintf(stderr, "  -f | --ref                Reference sequence used for mapping\n");
     fprintf(stderr, "                            (Only required for --dindel)\n");
     fprintf(stderr, "  -o | --out FILE           Output BAM file [- = stdout = default]\n");
     fprintf(stderr, "       --verbose            Be verbose\n");
     fprintf(stderr, "\n");
     fprintf(stderr,
             "The preferred way of inserting indel qualities should be via GATK's BQSR (>=2)" \
             " If that's not possible, use this subcommand.\n"  \
             "The command has two modes: 'uniform' and 'dindel':\n" \
             "- 'uniform' will assign a given value uniformly, whereas\n"  \
             "- 'dindel' will insert indel qualities based on Dindel (PMID 20980555).\n" \
             "Both will overwrite any existing values.\n");
     fprintf(stderr, "Do not realign your BAM file afterwards!\n");
     fprintf(stderr, "\n");
}


void idq_from_arg(int *iq, int *dq, const char *arg) 
{
     char *arg2 = strdup(arg);
     char *cpos = strchr(arg2, ',');
     if (cpos) {
          (*dq) = atoi(cpos+1);
          (*cpos) = '\0';
          (*iq) = atoi(arg2);
     } else {
          (*iq) = (*dq) = atoi(arg);
     }
     free(arg2);
}

int main_indelqual(int argc, char *argv[])
{
     char *bam_in = NULL;
     char *bam_out = NULL; /* - == stdout */
     char *ref = NULL;
     int c;
     static int dindel = 0;
     int uni_iq = -1;
     int uni_dq = -1;
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
               idq_from_arg(& uni_iq, & uni_dq, optarg);
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
          fprintf(stderr, "FATAL: Need exactly one BAM file as last argument\n");
          usage();
          return 1;
     }
     bam_in = (argv + optind + 1)[0];
     if ((0 != strcmp(bam_in, "-")) && ! file_exists(bam_in)) {
          LOG_FATAL("BAM file %s does not exist. Exiting...\n", bam_in);
          return -1;
     }

     if (! bam_out) {
          bam_out = malloc(2 * sizeof(char));
          strcpy(bam_out, "-");
     }

     LOG_DEBUG("uni_iq=%d\n", uni_iq);
     LOG_DEBUG("uni_dq=%d\n", uni_dq);
     LOG_DEBUG("bam_in=%s\n", bam_in);
     LOG_DEBUG("bam_out=%s\n", bam_out);
     LOG_DEBUG("ref=%s\n", ref);

     if ((uni_iq != -1 && uni_dq == -1)
         ||
         (uni_iq == -1 && uni_dq != -1)) {
          LOG_FATAL("internal logic error: uni_iq=%d uni_dq=%d\n", uni_iq, uni_dq);
          exit(1);
     }

     if (uni_iq != -1 && uni_dq != -1) {
          if (dindel) {
               LOG_FATAL("%s\n", "Can't insert both, uniform and dindel qualities");
               return -1;
          }
          return add_uniform(bam_in, bam_out, uni_iq, uni_dq);

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
