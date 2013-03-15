/* -*- c-file-style: "k&r" -*-
 *
 * This file is largely based on lofreq_snpcaller.c
 *
 * FIXME missing license
 *
 */
#include <math.h>
#include <stdio.h>
#include <unistd.h>
#include <ctype.h>
#include <string.h>
#include <errno.h>
#include <assert.h>
#include <stdarg.h>
#include <limits.h>
#include <float.h>
#include <getopt.h>
#include <stdlib.h>

#include "faidx.h"
#include "sam.h"
#include "kstring.h"

/* from bedidx.c */
void *bed_read(const char *fn);
void bed_destroy(void *_h);
int bed_overlap(const void *_h, const char *chr, int beg, int end);

#include "snpcaller.h"
#include "vcf.h"
#include "bam2depth.h"
#include "fet.h"
#include "utils.h"
#include "log.h"
#include "plp.h"



static void
usage(const mplp_conf_t *mplp_conf)
{
     fprintf(stderr, "Usage: %s call [options] in.bam\n\n", PACKAGE);
     fprintf(stderr, "Options:\n");
     /* generic */
     fprintf(stderr, "          --verbose           be verbose\n");
     fprintf(stderr, "          --debug             enable debugging\n");
     /* regions */
     fprintf(stderr, "       -r|--region STR        region in which pileup should be generated [null]\n");
     fprintf(stderr, "       -l|--bed FILE          list of positions (chr pos) or regions (BED) [null]\n");
     /*  */
     fprintf(stderr, "       -d|--maxdepth INT      max per-BAM depth to avoid excessive memory usage [%d]\n", mplp_conf->max_depth);
     fprintf(stderr, "       -f|--reffa FILE        faidx indexed reference sequence file [null]\n");
     /* */
     fprintf(stderr, "       -o|--out FILE          vcf output file [- = stdout]\n");
     /* base call quality and baq */
     fprintf(stderr, "       -q|--min-bq INT        skip any base with baseQ smaller than INT [%d]\n", mplp_conf->min_bq);
     fprintf(stderr, "       -B|--no-baq            disable BAQ computation\n");
     /* fprintf(stderr, "       -E           extended BAQ for higher sensitivity but lower specificity\n"); */
     /* mapping quality */
     fprintf(stderr, "       -m|--min_mq INT        skip alignments with mapQ smaller than INT [%d]\n", mplp_conf->min_mq);
     fprintf(stderr, "       -M|--max_mq INT        cap mapping quality at INT [%d]\n", mplp_conf->max_mq);
     fprintf(stderr, "       -J|--no-mq             don't merge mapQ into baseQ: P_e = P_mq + (1-P_mq) P_bq\n");
     /* stats */
     fprintf(stderr, "       -6|--illumina-1.3      assume the quality is Illumina-1.3-1.7/ASCII+64 encoded\n");
     fprintf(stderr, "       -A|--use-orphan        count anomalous read pairs\n");
}
/* usage() */



void
plp_summary(const plp_col_t *plp_col, void* conf) 
{
      FILE* stream = stdout;
     fprintf(stream, "%s\t%d\t%c\t%c\n",
             plp_col->target, plp_col->pos+1, plp_col->ref_base, plp_col->cons_base);
     
     LOG_FIXME("%s\n", "unfinished");
}

int 
main_plp_summary(int argc, char *argv[])
{
     /* based on bam_mpileup() */
     int c, i;
     static int use_orphan = 0;
     char *bam_file;
     char *bed_file = NULL;
     mplp_conf_t mplp_conf;
     void (*plp_proc_func)(const plp_col_t*, void *) = &plp_summary;

     LOG_FIXME("%s\n", "- Proper source qual use missing");
     LOG_FIXME("%s\n", "- Indel handling missing");
     LOG_FIXME("%s\n", "- Implement routine test against old SNV caller");
     LOG_FIXME("%s\n", "- Test actual SNV and SB values for both types of SNVs");

    memset(&mplp_conf, 0, sizeof(mplp_conf_t));
    /* default pileup options */
    mplp_conf.max_mq = 255; /* 60 */
    mplp_conf.min_bq = 3; /* 13 */
    mplp_conf.capQ_thres = 0;
    mplp_conf.max_depth = 1000000; /* 250 */
    mplp_conf.flag = MPLP_NO_ORPHAN | MPLP_REALN | MPLP_EXT_BAQ;


    /* FIXME getopt should be replaced with something more sensible like
     * argtable2 ot Gopt. Otherwise there's always the risk between
     * incosistent long opt, short opt and usage */
    while (1) {
         static struct option long_opts[] = {
              /* see usage sync */
              {"verbose", no_argument, &verbose, 1},
              {"debug", no_argument, &debug, 1},

              {"region", required_argument, NULL, 'r'},
              {"bed", required_argument, NULL, 'l'},
              
              {"maxdepth", required_argument, NULL, 'd'},
              {"reffa", required_argument, NULL, 'f'},
               
              {"out", required_argument, NULL, 'o'},

              {"min-bq", required_argument, NULL, 'q'},
              {"no-baq", no_argument, NULL, 'B'},
              /*{"ext-baq", required_argument, NULL, 'E'},*/
                   
              {"min-mq", required_argument, NULL, 'm'},
              {"max-mq", required_argument, NULL, 'M'},

              {"illumina-1.3", no_argument, NULL, 'I'},
              {"use-orphan", no_argument, &use_orphan, 1},

              {"help", no_argument, NULL, 'h'},

              {0, 0, 0, 0} /* sentinel */
         };

         /* see usage sync */
         static const char *long_opts_str = "r:l:d:f:o:q:Bm:M:Ih";

         /* getopt_long stores the option index here. */
         int long_opts_index = 0;
         c = getopt_long(argc, argv, long_opts_str, long_opts, & long_opts_index);
         if (c == -1) {
              break;
         }

         switch (c) {
         /* see usage sync */
         case 'r': mplp_conf.reg = strdup(optarg); break; /* FIXME you can enter lots of invalid stuff and libbam won't complain. add checks here? */
         case 'l': 
              mplp_conf.bed = bed_read(optarg); 
              bed_file = strdup(optarg);
              break;
              
         case 'd': mplp_conf.max_depth = atoi(optarg); break;
         case 'f':
              mplp_conf.fa = strdup(optarg);
              mplp_conf.fai = fai_load(optarg);
              if (mplp_conf.fai == 0) 
                   return 1;
              break;
         case 'q': mplp_conf.min_bq = atoi(optarg); break;
         case 'B': mplp_conf.flag &= ~MPLP_REALN; break;
         /* case 'E': mplp.flag |= MPLP_EXT_BAQ; break; */
              
         case 'm': mplp_conf.min_mq = atoi(optarg); break;
         case 'M': mplp_conf.max_mq = atoi(optarg); break;
         case 'I': mplp_conf.flag |= MPLP_ILLUMINA13; break;

         case 'h': usage(& mplp_conf); exit(0); /* WARN: not printing defaults if some args where parsed */
         case '?': LOG_FATAL("%s\n", "unrecognized arguments found. Exiting...\n"); exit(1);
#if 0
         case 0:  fprintf(stderr, "ERROR: long opt (%s) not mapping to short option. Exiting...\n", long_opts[long_opts_index].name); exit(1);
#endif
         default:
              break;
         }
    }
    if (use_orphan) {
         mplp_conf.flag &= ~MPLP_NO_ORPHAN;
    }
    mplp_conf.cmdline[0] = '\0';
    for (i=0; i<argc; i++) {
         strncat(mplp_conf.cmdline, argv[i], sizeof(mplp_conf.cmdline)-strlen(mplp_conf.cmdline)-2);
         strcat(mplp_conf.cmdline, " ");
    }

    if (argc == 1) {
        fprintf(stderr, "\n");
        usage(& mplp_conf);
        return 1;
    }
    if (1 != argc - optind) {
        fprintf(stderr, "Need exactly one BAM file as last argument\n");
        return(EXIT_FAILURE);
    }
    bam_file = (argv + optind)[0];


    if (debug) {
         dump_mplp_conf(& mplp_conf, stderr);
    }

    /* FIXME: implement logic_check_opts() */
    assert(mplp_conf.min_mq <= mplp_conf.max_mq);
    assert(! (mplp_conf.bed && mplp_conf.reg));
   
    (void) mpileup(&mplp_conf, (void*)plp_proc_func, NULL,
                   1, (const char **) argv + optind);


    free(mplp_conf.reg); 
    free(mplp_conf.fa);
    if (mplp_conf.fai) {
         fai_destroy(mplp_conf.fai);
    }
    free(bed_file);
    if (mplp_conf.bed) {
         bed_destroy(mplp_conf.bed);
    }
    LOG_VERBOSE("%s\n", "Successful exit.");
    return 0;
}
/* main_call */


#ifdef MAIN_TEST


int main()
{     
     fprintf(stderr, "WARNING: main replaced with test function. just monkeying around\n");

     if (1) {
          char test_nucs[] = "ACGTNRYacgtnryZ\0";
          int i;
          
          for (i=0; i<strlen(test_nucs); i++) {
               printf("%d %c - %d - %d - %d\n", i, test_nucs[i],
                      bam_nt16_table[(int)test_nucs[i]],
                      bam_nt16_nt4_table[bam_nt16_table[(int)test_nucs[i]]],
                      bam_nt4_table[(int)test_nucs[i]]);
          }
     }

     if (1) {
          int quals[] = {30, 30, 30, 30, 30, 30, 30, 30, 30, 30,
                         30, 30, 30, 30, 30, 30, 30, 30, 30, 30,
                         30, 30, 30, 30, 30, 30, 30, 30, 30, 30,
                         30, 30, 30, 30, 30, 30, 30, 30, 30, 30,
                         30, 30, 30, 30, 30, 30, 30, 30, 30, 30};
          int n_quals = 50;
          int n_mismatches = 1;
          double *probvec;
          double src_pvalue;
          int src_qual; 
          int i;

          qsort(quals, n_quals, sizeof(int), int_cmp);

          for (n_mismatches=0; n_mismatches<n_quals/2; n_mismatches++) {
               probvec = poissbin(&src_pvalue, quals,
                                  n_quals, n_mismatches, 1.0, 0.05);
               
               if (src_pvalue>1.0) {/* DBL_MAX is default return value */
                    src_pvalue = 1.0;/*-DBL_EPSILON;*/
               }
               /* src_pvalue: what's the chance of seeing n_mismatches or more
                * given quals? or: how likely is this read from the genome.
                * 1-src_value = prob read is not from genome
                */
               LOG_FIXME("Orig src_pv = %f", src_pvalue);
               src_pvalue = 1.0-src_pvalue;
               free(probvec);
               
               src_qual =  PROB_TO_PHREDQUAL(src_pvalue);
               fprintf(stderr, "| src_pv = %f = Q%d for %d/%d mismatches. All quals: ", 
                       src_pvalue, src_qual, n_mismatches, n_quals);
               for (i=0; i<n_quals; i++) {
                    fprintf(stderr, " %d", quals[i]);
               }
               fprintf(stderr, "\n");
          }
     }
     return 0;
}

#endif

