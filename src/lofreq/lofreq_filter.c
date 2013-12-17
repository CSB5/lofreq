/* -*- c-file-style: "k&r"; indent-tabs-mode: nil; -*- */


/*********************************************************************
 *
 * FIXME update license
 *
 *********************************************************************/


#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <assert.h>
#include <stdarg.h>
#include <getopt.h>
#include <stdlib.h>

/* lofreq includes */
#include "lofreq_filter.h"
#include "vcf.h"
#include "log.h"
#include "utils.h"
#include "multtest.h"



#if 1
#define MYNAME "lofreq filter"
#else
#define MYNAME PACKAGE
#endif

#define DEFAULT_ALPHA 0.05


typedef struct {
     int min;
     int max;
} cov_filter_t;

typedef struct {
     float min;
     float max;
} af_filter_t;

typedef struct {
     int thresh;/* use if > 0; otherwise use multiple testing correction that's if >0 */
     int mtc_type;/* holm; holmbonf; fdr; none */
     double alpha;
     int ntests;
} sb_filter_t;

typedef struct {
     int thresh;/* use if > 0; otherwise use multiple testing correction that's if >0 */
     int mtc_type;/* holm; holmbonf; fdr; none */
     double alpha;
     int ntests;
} snvqual_filter_t;

typedef struct {
     vcf_file_t vcf_in;
     vcf_file_t vcf_out;
     int print_only_passed;

     /* each allowed to be NULL if not set */
     cov_filter_t cov_filter;
     af_filter_t af_filter;
     sb_filter_t sb_filter;
     snvqual_filter_t snvqual_filter;
} filter_conf_t;


void
dump_filter_conf(filter_conf_t *cfg)
{
     fprintf(stderr, "filter_conf:\n");
     fprintf(stderr, "  print_only_passed=%d\n", cfg->print_only_passed);

     fprintf(stderr, "  cov_filter min=%d max=%d\n",
             cfg->cov_filter.min, cfg->cov_filter.max);
     fprintf(stderr, "  af_filter min=%f max=%f\n", 
             cfg->af_filter.min, cfg->af_filter.max);
     fprintf(stderr, "  sb_filter thresh=%d mtc_type=%d|%s alpha=%f ntests=%d\n", 
             cfg->sb_filter.thresh, cfg->sb_filter.mtc_type, mtc_type_str[cfg->sb_filter.mtc_type], cfg->sb_filter.alpha, cfg->sb_filter.ntests);
     fprintf(stderr, "  snvqual_filter thresh=%d mtc_type=%d|%s alpha=%f ntests=%d\n", 
             cfg->snvqual_filter.thresh, cfg->snvqual_filter.mtc_type, mtc_type_str[cfg->snvqual_filter.mtc_type], cfg->snvqual_filter.alpha, cfg->snvqual_filter.ntests);
}


static void
usage(const filter_conf_t* filter_conf)
{
     fprintf(stderr, "%s: Filter variant parsed from vcf file\n\n", MYNAME);
     fprintf(stderr, "Usage: %s [options] -i input.vcf -o output.vcf\n", MYNAME);

     fprintf(stderr,"Options:\n");
     fprintf(stderr, "  Files:\n"); 
     fprintf(stderr, "  -i | --in FILE              VCF input file (gzip supported)\n");
     fprintf(stderr, "  -o | --out FILE             VCF output file (default: - for stdout; gzip supported).\n");

     fprintf(stderr, "  Coverage:\n");
     fprintf(stderr, "  -v | --cov-min INT          Minimum coverage allowed (<1=off)\n");
     fprintf(stderr, "  -V | --cov-max INT          Maximum coverage allowed (<1=off)\n");

     fprintf(stderr, "  Allele Frequency (neg. values = off):\n");
     fprintf(stderr, "  -a | --af-min FLOAT         Maximum allele freq allowed (<1=off)\n");
     fprintf(stderr, "  -A | --af-max FLOAT         Minimum allele freq allowed (<1=off)\n");

     fprintf(stderr, "  Strand Bias:\n");
     fprintf(stderr, "  -B | --sb-thresh INT        Maximum phred-value allowed. Conflicts with -b.\n");
     fprintf(stderr, "  -b | --sb-mtc STRING        Multiple testing correction type. One of bonf, holm or fdr. Conflicts with -B\n");
     fprintf(stderr, "  -c | --sb-alpha FLOAT       Multiple testing correcion pvalue threshold\n");

     fprintf(stderr, "  SNV quality:\n");
     fprintf(stderr, "  -Q  | --snvqual-thresh INT  Maximum phred-value allowed. Conflicts with -q\n");
     fprintf(stderr, "  -q  | --snvqual-mtc STRING  Multiple testing correction type. One of bonf, holm or fdr. Conflicts with -Q\n");
     fprintf(stderr, "  -r  | --snvqual-alpha FLOAT Multiple testing correcion pvalue threshold\n");
     fprintf(stderr, "  -s  | --snvqual-ntests INT  Multiple testing correcion pvalue threshold\n");

     fprintf(stderr, "  Misc.:\n");
     fprintf(stderr, "       --verbose              Be verbose\n");
     fprintf(stderr, "       --debug                Enable debugging\n");
     fprintf(stderr, "       --only-passed          Only output passed variants\n");
}
/* usage() */


int 
main_filter(int argc, char *argv[])
{
     filter_conf_t cfg;
     char *vcf_in = NULL, *vcf_out = NULL;
     static int print_only_passed = 0;
     char *vcf_header = NULL;
     var_t **vars = NULL;
     long int num_vars = 0;
     long int vars_size = 0;
     long int i;

     /* default filter options */
     memset(&cfg, 0, sizeof(filter_conf_t));
     cfg.cov_filter.min = cfg.cov_filter.max = -1;
     cfg.af_filter.min = cfg.af_filter.max = -1;
     cfg.sb_filter.alpha = DEFAULT_ALPHA;
     cfg.snvqual_filter.alpha = DEFAULT_ALPHA;


    /* keep in sync with long_opts_str and usage 
     *
     * getopt is a pain in the whole when it comes to syncing of long
     * and short args and usage. check out gopt, libcfu...
     */
    while (1) {
         int c;
         static struct option long_opts[] = {
              /* see usage sync */
              {"verbose", no_argument, &verbose, 1},
              {"debug", no_argument, &debug, 1},
              {"only-passed", no_argument, &print_only_passed, 1},

              {"help", no_argument, NULL, 'h'},
              {"in", required_argument, NULL, 'i'},
              {"out", required_argument, NULL, 'o'},

              {"cov-min", required_argument, NULL, 'v'},
              {"cov-max", required_argument, NULL, 'V'},

              {"af-min", required_argument, NULL, 'a'},
              {"af-max", required_argument, NULL, 'A'},

              {"sb-thresh", required_argument, NULL, 'B'},
              {"sb-mtc", required_argument, NULL, 'b'},
              {"sb-alpha", required_argument, NULL, 'c'},

              {"snvqual-thresh", required_argument, NULL, 'Q'},
              {"snvqual-mtc", required_argument, NULL, 'q'},
              {"snvqual-alpha", required_argument, NULL, 'r'},
              {"snvqual-ntests", required_argument, NULL, 's'},

              {0, 0, 0, 0} /* sentinel */
         };

         /* keep in sync with long_opts and usage */
         static const char *long_opts_str = "hi:o:v:V:a:A:B:b:c:Q:q:r:s:"; 

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
              usage(& cfg); 
              return 0;

         case 'i': 
              vcf_in = strdup(optarg);
              break;
         case 'o':
              if (0 != strcmp(optarg, "-")) {
                   if (file_exists(optarg)) {
                        LOG_FATAL("Cowardly refusing to overwrite file '%s'. Exiting...\n", optarg);
                        return 1;
                   }
              }
              vcf_out = strdup(optarg);
              break;

         case 'v': 
              cfg.cov_filter.min = atoi(optarg);
              break;
         case 'V': 
              cfg.cov_filter.max = atoi(optarg);
              break;

         case 'a': 
              cfg.af_filter.min = strtof(optarg, NULL);
              break;
         case 'A': 
              cfg.af_filter.max = strtof(optarg, NULL);
              break;

         case 'B':
              cfg.sb_filter.thresh = atoi(optarg);
              break;
         case 'b':
              cfg.sb_filter.mtc_type = mtc_str_to_type(optarg);
              if (-1 == cfg.sb_filter.mtc_type) {
                   LOG_FATAL("Unknown multiple testing correction type '%s' for strandbias filtering\n", optarg);
                   return -1;
              }
              break;
         case 'c':
              cfg.sb_filter.alpha = strtof(optarg, NULL);
              break;

         case 'Q':
              cfg.snvqual_filter.thresh = atoi(optarg);
              break;
         case 'q':
              cfg.snvqual_filter.mtc_type = mtc_str_to_type(optarg);
              if (-1 == cfg.snvqual_filter.mtc_type) {
                   LOG_FATAL("Unknown multiple testing correction type '%s' for snv quality filtering\n", optarg);
                   return -1;
              }
              break;
         case 'r':
              cfg.snvqual_filter.alpha = strtof(optarg, NULL);
              break;
         case 's':
              cfg.snvqual_filter.ntests = atoi(optarg);
              break;

         case '?': 
              LOG_FATAL("%s\n", "Unrecognized argument found. Exiting...\n"); 
              return 1;

         default:
              break;
         }
    }
    cfg.print_only_passed = print_only_passed;

    if (0 != argc - optind - 1) {/* FIXME needed at all? */
         LOG_FATAL("%s\n", "Unrecognized argument found. Exiting...\n");
         return 1;
    }

    /* logic check of command line parameters
     */
    if (cfg.cov_filter.max > 0 &&  cfg.cov_filter.max < cfg.cov_filter.min) {
         LOG_FATAL("%s\n", "Invalid coverage-filter settings");
         return 1;
    }
    if ((cfg.af_filter.max > 0 && cfg.af_filter.max < cfg.af_filter.min) ||
        (cfg.af_filter.max > 1.0)) {
         LOG_FATAL("%s\n", "Invalid AF-filter settings");
         return 1;
    }

    if (cfg.sb_filter.thresh && cfg.sb_filter.mtc_type != MTC_NONE) {
         LOG_FATAL("%s\n", "Can't use fixed strand-bias threshold *and* multiple testing correction.");
         return 1;
    }
    if (cfg.snvqual_filter.thresh && cfg.snvqual_filter.mtc_type != MTC_NONE) {
         LOG_FATAL("%s\n", "Can't use fixed SNV quality threshold *and* multiple testing correction.");
         return 1;
    }

    LOG_FIXME("add more checks esp. for sb and snvqual. Also init values\n");

    if (argc == 2) {
        fprintf(stderr, "\n");
        usage(& cfg);
        return 1;
    }

    dump_filter_conf(& cfg);

    /* missing file args default to stdin and stdout 
     */
    if  (! vcf_in) {
         vcf_in = malloc(2 * sizeof(char));
         strcpy(vcf_in, "-");
    }
    if  (! vcf_out) {
         vcf_out = malloc(2 * sizeof(char));
         strcpy(vcf_out, "-");
    }
    LOG_DEBUG("vcf_in=%s vcf_out=%s\n", vcf_in, vcf_out);

    /* open vcf files
     */
    if (vcf_file_open(& cfg.vcf_in, vcf_in, 
                      HAS_GZIP_EXT(vcf_in), 'r')) {
         LOG_ERROR("Couldn't open %s\n", vcf_in);
         return 1;
    }
    if (vcf_file_open(& cfg.vcf_out, vcf_out, 
                      HAS_GZIP_EXT(vcf_out), 'w')) {
         LOG_ERROR("Couldn't open %s\n", vcf_out);
         return 1;
    }
    free(vcf_in);
    free(vcf_out);

    /* print header
     */
    if (0 !=  vcf_parse_header(&vcf_header, & cfg.vcf_in)) {
         LOG_WARN("%s\n", "vcf_parse_header() failed");
         if (vcf_file_seek(& cfg.vcf_in, 0, SEEK_SET)) {
              LOG_FATAL("%s\n", "Couldn't rewind file to parse variants"
                        " after header parsing failed");
              return -1;
         }
    } 

    LOG_FIXME("%s\n", "Add filter to header");
    vcf_write_header(& cfg.vcf_out, vcf_header);
    free(vcf_header);

    /* read variants. since many filters perform multiple testing
     * correction and therefore need to look at all variants we keep
     * it simple and load them all into memory.
     */
    num_vars = 0;
    while (1) {
         var_t *var;
         int rc;
         vcf_new_var(&var);
         rc = vcf_parse_var(& cfg.vcf_in, var);
         if (-1 == rc) {
              LOG_FATAL("%s\n", "Error while parsing vcf-file");
              exit(1);
         }
                   
         if (1 == rc) {/* EOF */
              free(var);
              break;
         }

         /* read all in no matter if already filtered. we keep adding to filters */         
         num_vars +=1;
         if (num_vars >= vars_size) {
              const long incr = 128;
              vars = realloc(vars, (vars_size+incr) * sizeof(var_t**));
              vars_size += incr;
         }
         vars[num_vars-1] = var;
#ifdef TRACE
         {
              char *key;
              vcf_var_key(&key,  vars[num_vars-1]);
              fprintf(stderr, "storing var %ld+1: %s\n", num_vars, key);
              free(key);
         }
#endif
    }
    if (num_vars) {
         vars = realloc(vars, (num_vars * sizeof(var_t**)));
    }
    vcf_file_close(& cfg.vcf_in);
    LOG_FIXME("%s\n", "init ntests in sb_ and snvqual_ if type != MTC_NONE");

    LOG_FIXME("Got %d vars\n", num_vars);
    for (i=0; i<num_vars; i++) {
         var_t *v = vars[i];
         if (cfg.print_only_passed && ! (VCF_VAR_PASSES(v))) {
              continue;
         }
         vcf_write_var(& cfg.vcf_out, v);
    }

    LOG_FIXME("%s\n", "Implement lazy option:\n" \
              " - suck all into mem" \
              " - apply filter one by one\n" \
              " - final iteration: set pass to any without filter tag");
    LOG_FIXME("%s\n", "unfinished");

    vcf_file_close(& cfg.vcf_out);

    LOG_VERBOSE("%s\n", "Successful exit.");
    return 0;
}
/* main_filter */


/* cc lofreq_filter.c -o lofreq_filter -I../lofreq_core -I../uthash/ ../lofreq_core/liblofreq_core.a   -lz -DMAIN_FILTER */
#ifdef MAIN_FILTER
int 
main(int argc, char *argv[])
{
     return main_filter(argc+1, argv-1);
}
#endif
