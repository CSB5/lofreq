/* -*- c-file-style: "k&r"; indent-tabs-mode: nil; -*- */


/*********************************************************************
 *
 * FIXME update license
 *
 *********************************************************************/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>

/* samtools includes */
#include "sam.h"
#include "faidx.h"
/* from bedidx.c */
void *bed_read(const char *fn);
void bed_destroy(void *_h);
int bed_overlap(const void *_h, const char *chr, int beg, int end);

/* lofreq includes */
#include "log.h"
#include "utils.h"
#include "samutils.h"


#if 1
#define MYNAME "lofreq bam-stats"
#else
#define MYNAME PACKAGE
#endif



typedef struct {
     int min_mq;
     int min_bq;
     char *fa;
     faidx_t *fai;
     void *bed;
     int samflags_on;
     int samflags_off;
} bamstats_conf_t;

/* sam_view.c:__g_skip_aln */
static inline int skip_aln(const bam_header_t *h, const bam1_t *b,
                           const int min_mq, const int flag_on, const int flag_off, void *bed)
{
     if (b->core.qual < min_mq || ((b->core.flag & flag_on) != flag_on) || (b->core.flag & flag_off)) {
          /*fprintf(stderr, "Skipping because of flag or min_mq\n");*/
          return 1;
     }
     if (bed && b->core.tid >= 0 && !bed_overlap(bed, h->target_name[b->core.tid], b->core.pos, bam_calend(&b->core, bam1_cigar(b)))) {
          /*fprintf(stderr, "Skipping because of bed: h->target_name[b->core.tid=%d] = %s; b->core.pos = %d\n", b->core.tid, h->target_name[b->core.tid], b->core.pos);*/
          return 2;
     }
     /*fprintf(stderr, "Not skipping\n");*/
     return 0;
}



static void
usage(bamstats_conf_t *bamstats_conf)
{
     fprintf(stderr, "%s: Compiles statistics from BAM files\n\n", MYNAME);

     fprintf(stderr,"Usage: %s [options] in.bam\n\n", MYNAME);
     fprintf(stderr,"Options:\n");
     fprintf(stderr, "       --verbose       Be verbose\n");
     fprintf(stderr, "       --debug         Enable debugging\n");
     fprintf(stderr, "  -l | --bed FILE      List of positions (chr pos) or regions (BED) [null]\n");
     fprintf(stderr, "  -f | --reffa FILE    Indexed reference fasta file (gzip supported) [null]\n");
     fprintf(stderr, "  -o | --out FILE      Write stats to this output file [- = stdout]\n");
     fprintf(stderr, "  -q | --min-bq INT    Ignore any base with baseQ smaller than INT [%d]\n", bamstats_conf->min_bq);
     fprintf(stderr, "  -m | --min-mq INT    Ignore reads with mapQ smaller than INT [%d]\n", bamstats_conf->min_mq);
}
/* usage() */




int 
main_bamstats(int argc, char *argv[])
{
     char *outfile = NULL;
     char *bamfile = NULL;
     char *bedfile = NULL;
     samfile_t *sam =  NULL;
     long int num_good_reads = 0;
     long int num_ign_reads = 0;
     int rc = 0;
     bamstats_conf_t bamstats_conf;

     memset(&bamstats_conf, 0, sizeof(bamstats_conf_t));
     bamstats_conf.min_mq = 13;
     bamstats_conf.min_bq = 3;
     /* will skip read if any of the following is set */
     bamstats_conf.samflags_off = 0;
     bamstats_conf.samflags_off |= 0x4; /* segment unmapped */
     bamstats_conf.samflags_off |= 0x100; /* secondary alignment */
     bamstats_conf.samflags_off |= 0x200; /* not passing quality controls */
     bamstats_conf.samflags_off |= 0x400; /* PCR or optical duplicate */
     bamstats_conf.samflags_off |= 0x800; /* supplementary alignment */

     /* FIXME enable BAQ on request ? */

     /* keep in sync with long_opts_str and usage 
     *
     * getopt is a pain in the whole when it comes to syncing of long
     * and short args and usage. check out gopt, libcfu
     */
     while (1) {
          int c;
          static struct option long_opts[] = {
               /* see usage sync */               
               {"bed", required_argument, NULL, 'l'},
               {"reffa", required_argument, NULL, 'f'},
               {"out", required_argument, NULL, 'o'},
               {"min-bq", required_argument, NULL, 'q'},
               {"min-mq", required_argument, NULL, 'm'},

               {"help", no_argument, NULL, 'h'},
               {"verbose", no_argument, &verbose, 1},
               {"debug", no_argument, &debug, 1},

               {0, 0, 0, 0} /* sentinel */
          };
          
          /* keep in sync with long_opts and usage */
          static const char *long_opts_str = "hl:f:o:q:m:"; 
          
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
               usage(&bamstats_conf); 
               return 0;

          case 'l': 
              bedfile = strdup(optarg);
              break;

          case 'f':
               bamstats_conf.fa = strdup(optarg);
               bamstats_conf.fai = fai_load(optarg);
               if (bamstats_conf.fai == 0)  {
                    free(bamstats_conf.fa);
                    return 1;
               }
               break;

          case 'o':
               if (0 != strcmp(optarg, "-")) {
                    if (file_exists(optarg)) {
                         LOG_FATAL("Cowardly refusing to overwrite file '%s'. Exiting...\n", optarg);
                         return 1;
                    }
               } 
               outfile = strdup(optarg);
               break;
               
         case 'q': 
              bamstats_conf.min_bq = atoi(optarg); 
              break;

          case 'm': 
              bamstats_conf.min_mq = atoi(optarg); 
              break;

         case '?': 
               LOG_FATAL("%s\n", "Unrecognized arguments found. Exiting...\n"); 
               return 1;
               
          default:
               break;
          }
     }
     
     
     if (argc == 2) {
          fprintf(stderr, "\n");
          usage(&bamstats_conf);
          return 1;
     }
     
     if (1 != argc - optind - 1) {
          LOG_FATAL("%s\n\n", "Need exactly one BAM file as last argument");
          return 1;
     }
     bamfile = (argv + optind + 1)[0];
     if (!file_exists(bamfile)) {
          LOG_FATAL("BAM file %s does not exist.\n\n", bamfile);
          return -1;
     }
     
     if (NULL == bamstats_conf.fa) {
          LOG_FATAL("%s\n\n", "ERROR: Missing reference fasta argument");
          usage(&bamstats_conf);
          return 1;
     }

     if (bedfile) {
          LOG_VERBOSE("%s\n", "Bed routines don't make use of indexing, i.e. read BAM file fully and are therefore slow."); /* FIXME */
          bamstats_conf.bed = bed_read(bedfile);
     }


    
     if ((sam = samopen(bamfile, "rb", NULL)) == 0) {
          LOG_FATAL("Failed to open \"%s\" for reading.\n", bamfile);
          rc = 1;

     } else {
          /* partially based on 0.1.18 sam_view.c */
          bam1_t *b = bam_init1();
          int r;
          char *target_name = NULL;
          char *ref = NULL;

          while ((r = samread(sam, b)) >= 0) { // read one alignment from `in'
               int counts[MAX_COUNT_IDX];
               int ref_len = -1;

               if (skip_aln(sam->header, b, bamstats_conf.min_mq, 
                              bamstats_conf.samflags_on, bamstats_conf.samflags_off,
                              bamstats_conf.bed)) {
                    num_ign_reads += 1;
                    continue;
               }

               /* reload ref if necessary */
               if (target_name == NULL || strcmp(target_name, sam->header->target_name[b->core.tid]) != 0) {
                    if (ref) {
                         free(ref);
                    }
                    LOG_FIXME("%s\n", "reloading ref");
                    target_name = sam->header->target_name[b->core.tid];
                    ref = faidx_fetch_seq(bamstats_conf.fai, target_name,
                                          0, 0x7fffffff, &ref_len);
               }
               if (count_matches(counts, b, ref, bamstats_conf.min_mq)) {
                    LOG_WARN("%s\n", "count_matches failed on read. ignoring"); /* FIXME print read */
                    continue;
               }

               num_good_reads += 1;

               LOG_FIXME("good/ign=%u/%u: m=%d mm=%d i=%d d=%d\n", num_good_reads, num_ign_reads,
                         counts[MATCH_COUNT_IDX], counts[MISMATCH_COUNT_IDX], counts[INS_COUNT_IDX], counts[DEL_COUNT_IDX]);
          }
          if (r < -1) {
               LOG_FATAL("%s\n", "BAM file is truncated.\n");
               rc = 1;
          }
          bam_destroy1(b);
     }


     /* cleanup */

     samclose(sam);

     free(outfile);

     free(bamstats_conf.fa);
     if (bamstats_conf.fai) {
         fai_destroy(bamstats_conf.fai);
     }

     free(bedfile);
     if (bamstats_conf.bed) {
          bed_destroy(bamstats_conf.bed);
     }

     if (0==rc) {
          LOG_VERBOSE("%s\n", "Successful exit.");
     }

     return rc;
}
/* main_bamstats */
