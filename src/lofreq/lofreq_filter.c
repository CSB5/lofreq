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



#if 1
#define MYNAME "lofreq filter"
#else
#define MYNAME PACKAGE
#endif


typedef struct {
     vcf_file_t vcf_in;
     vcf_file_t vcf_out;
     int print_only_passed;
} filter_conf_t;




static void
usage(const filter_conf_t* filter_conf)
{
     /* FIXME all */
     fprintf(stderr, "%s: Filter variant parsed from vcf file\n\n", MYNAME);
     fprintf(stderr, "Usage: %s [options] -i input.vcf -o output.vcf\n", MYNAME);

     fprintf(stderr,"Options:\n");
     fprintf(stderr, "       --verbose        Be verbose\n");
     fprintf(stderr, "       --debug          Enable debugging\n");
     fprintf(stderr, "       --only-passed    Only output passed variants\n");
     fprintf(stderr, "  -i | --in FILE        VCF input file (gzip supported)\n");
     fprintf(stderr, "  -o | --out FILE       VCF output file (default: - for stdout; gzip supported).\n");
}
/* usage() */




int 
main_filter(int argc, char *argv[])
{
     filter_conf_t filter_conf;
     char *vcf_in, *vcf_out;
     static int print_only_passed = 0;
     char *vcf_header = NULL;

     vcf_in = vcf_out = NULL;

     /* default filter options */
     memset(&filter_conf, 0, sizeof(filter_conf_t));


    /* keep in sync with long_opts_str and usage 
     *
     * getopt is a pain in the whole when it comes to syncing of long
     * and short args and usage. check out libcfu (also has hash
     * functions etc)
     */
    while (1) {
         int c;
         static struct option long_opts[] = {
              /* see usage sync */
              {"help", no_argument, NULL, 'h'},
              {"verbose", no_argument, &verbose, 1},
              {"debug", no_argument, &debug, 1},
              {"only-passed", no_argument, &print_only_passed, 1},

              {"in", required_argument, NULL, 'i'},
              {"out", required_argument, NULL, 'o'},

              {0, 0, 0, 0} /* sentinel */
         };

         /* keep in sync with long_opts and usage */
         static const char *long_opts_str = "hi:o:"; 

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
              usage(& filter_conf); 
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

         case '?': 
              LOG_FATAL("%s\n", "unrecognized arguments found. Exiting...\n"); 
              return 1;
         default:
              break;
         }
    }
    filter_conf.print_only_passed = print_only_passed;


    if (argc == 2) {
        fprintf(stderr, "\n");
        usage(& filter_conf);
        return 1;
    }


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

    /* open vcf files
     */
    if (vcf_file_open(& filter_conf.vcf_in, vcf_in, 
                      HAS_GZIP_EXT(vcf_in), 'r')) {
         LOG_ERROR("Couldn't open %s\n", vcf_in);
         return 1;
    }
    if (vcf_file_open(& filter_conf.vcf_out, vcf_out, 
                      HAS_GZIP_EXT(vcf_out), 'w')) {
         LOG_ERROR("Couldn't open %s\n", vcf_out);
         return 1;
    }

    free(vcf_in);
    free(vcf_out);


    /* print header
     */
    if (0 !=  vcf_parse_header(&vcf_header, & filter_conf.vcf_in)) {
         LOG_WARN("%s\n", "vcf_parse_header() failed");
         if (vcf_file_seek(& filter_conf.vcf_in, 0, SEEK_SET)) {
              LOG_FATAL("%s\n", "Couldn't rewind file to parse variants"
                        " after header parsing failed");
              return -1;
         }
    } else {
         vcf_write_header(& filter_conf.vcf_out, vcf_header);
         free(vcf_header);
    }



    /* read variants
     */
    while (1) {
         var_t *var;
         int rc;
         vcf_new_var(&var);
         rc = vcf_parse_var(& filter_conf.vcf_in, var);
         if (-1 == rc) {
              LOG_FATAL("%s\n", "Error while parsing vcf-file");
              exit(1);
         }
         if (1 == rc) {/* EOF */
              free(var);
              break;
         }
    }
    vcf_file_close(& filter_conf.vcf_in);
    vcf_file_close(& filter_conf.vcf_out);


    LOG_FIXME("%s\n", "PLEASE IMPLEMENT ME"); exit(1);


    LOG_VERBOSE("%s\n", "Successful exit.");
    return 0;
}
/* main_filter */

