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
#include "lofreq_vcfset.h"
#include "vcf.h"
#include "log.h"
#include "utils.h"



#if 1
#define MYNAME "lofreq vcfset"
#else
#define MYNAME PACKAGE
#endif


typedef enum {
     SETOP_UNKNOWN,
     SETOP_INTERSECT,
     SETOP_COMPLEMENT, 
     SETOP_UNION
} vcfset_op_t;

typedef struct {
     FILE *vcf_in1;
     FILE *vcf_in2;
     FILE *vcf_out;
     vcfset_op_t vcf_setop;
} vcfset_conf_t;




static void
usage(const vcfset_conf_t* vcfset_conf)
{
     fprintf(stderr, "%s: Perform set operations on two vcf files\n\n", MYNAME);
     fprintf(stderr, "Usage: %s [options] -a op -1 1.vcf -2 2.vcf \n", MYNAME);

     fprintf(stderr,"Options:\n");
     fprintf(stderr, "       --verbose        Be verbose\n");
     fprintf(stderr, "       --debug          Enable debugging\n");
     fprintf(stderr, "  -1 | --vcf1 FILE      1st VCF input file (gzip supported)\n");
     fprintf(stderr, "  -2 | --vcf2 FILE      2nd VCF input file (gzip supported)\n");
     fprintf(stderr, "  -o | --vcfout         VCF output file (- for stdout, which is default. gzip supported).\n"
             "                        Meta-data will be copied from vcf1\n");
     fprintf(stderr, "  -a | --action         Set operation to perform:\n"
             "                        intersect, complement or union.\n"
             "                        intersect = vcf1 AND vcf2.\n"
             "                        union = vcf1 OR vcf2.\n"
             "                        complement = vcf1 \\ vcf2.\n");
}
/* usage() */



int 
main_vcfset(int argc, char *argv[])
{
     vcfset_conf_t vcfset_conf;
     char *vcf_header = NULL;
     int rc = 0;
     int c, i;
     var_t **vars_vcf2 = NULL;
     long int num_vars_vcf1, num_vars_vcf2;

     num_vars_vcf1 = num_vars_vcf2 = 0;


     /* default uniq options */
     memset(&vcfset_conf, 0, sizeof(vcfset_conf_t));
     /* vcfset_conf.vcf_in1 = NULL; */
     /* vcfset_conf.vcf_in2 = NULL; */
     vcfset_conf.vcf_out = stdout;


    /* keep in sync with long_opts_str and usage 
     *
     * getopt is a pain in the whole when it comes to syncing of long
     * and short args and usage. check out libcfu (also has hash
     * functions etc)
     */
    while (1) {
         static struct option long_opts[] = {
              /* see usage sync */
              {"help", no_argument, NULL, 'h'},
              {"verbose", no_argument, &verbose, 1},
              {"debug", no_argument, &debug, 1},

              {"vcf1", required_argument, NULL, '1'},
              {"vcf2", required_argument, NULL, '2'},
              {"vcfout", required_argument, NULL, 'o'},
              {"action", required_argument, NULL, 'a'},

              {0, 0, 0, 0} /* sentinel */
         };

         /* keep in sync with long_opts and usage */
         static const char *long_opts_str = "h1:2:o:a:"; 

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
              usage(& vcfset_conf); 
              return 0;

         case '1': 
              if (! file_exists(optarg)) {
                   LOG_FATAL("Input file '%s' does not exist. Exiting...\n", optarg);
                   return 1;
              }
              vcfset_conf.vcf_in1 = fopen(optarg, "r");
              break;

         case '2': 
              if (! file_exists(optarg)) {
                   LOG_FATAL("Input file '%s' does not exist. Exiting...\n", optarg);
                   return 1;
              }
              vcfset_conf.vcf_in2 = fopen(optarg, "r");
              break;

         case 'o':
              if (0 != strcmp(optarg, "-")) {
                   if (file_exists(optarg)) {
                        LOG_FATAL("Cowardly refusing to overwrite file '%s'. Exiting...\n", optarg);
                        return 1;
                   }
                   vcfset_conf.vcf_out = fopen(optarg, "w");
              } else {
                   vcfset_conf.vcf_out = stdout;
                   
              }
              break;

         case 'a': 
              if (0 == strcmp(optarg, "intersect")) {
                   vcfset_conf.vcf_setop = SETOP_INTERSECT;
              } else if (0 == strcmp(optarg, "complement")) {
                   vcfset_conf.vcf_setop = SETOP_COMPLEMENT;
              } else if (0 == strcmp(optarg, "union")) {
                   vcfset_conf.vcf_setop = SETOP_UNION;
              } else {
                   LOG_FATAL("Unknown action '%s'. Exiting...\n", optarg);
                   return 1;
              }
              break;

         case '?': 
              LOG_FATAL("%s\n", "unrecognized arguments found. Exiting...\n"); 
              return 1;
         default:
              break;
         }
    }

    if (argc == 2) {
        fprintf(stderr, "\n");
        usage(& vcfset_conf);
        return 1;
    }


    if  (vcfset_conf.vcf_in1 == 0 ||
         vcfset_conf.vcf_in2 == 0 ||
         vcfset_conf.vcf_out == 0 ||
         vcfset_conf.vcf_setop == SETOP_UNKNOWN) {
         LOG_FATAL("Missing argument\n\n");
         usage(& vcfset_conf);
         return 1;
    }






    /* recipe: read B into memory and parse from A one by one
     * ======================================================
     */

    /* use meta-data/header of vcf_in1 for output
     */
    if (0 !=  vcf_parse_header(&vcf_header, vcfset_conf.vcf_in1)) {
         LOG_WARN("%s\n", "vcf_parse_header() failed");
         if (fseek(vcfset_conf.vcf_in1, 0, SEEK_SET)) {
              LOG_FATAL("%s\n", "Couldn't rewind file to parse variants"
                        " after header parsing failed");
              return -1;
         }
    } else {
         /* vcf_write_header would write *default* header */
         fprintf(vcfset_conf.vcf_out, "%s", vcf_header);
         free(vcf_header);
    }


    /* skip meta-data/header in vcf_in2
     */
    if (0 != vcf_skip_header(vcfset_conf.vcf_in2)) {
         LOG_FATAL("%s\n", "Failed to skip header in 2nd vcf files");
         return -1;
    }

    
    LOG_ERROR("parse from vcf2 one by one and save as hash\n"); /*  */

/*  
  see lofreq_uniq.c for vcf parsing examples
  
    hash = "%s %d %s %s" % (var.CHROM, var.POS, 
                            var.REF, ''.join(var.ALT))

    if (-1 == (num_vars_vcf2 = vcf_parse_vars(vcfset_conf.vcf_in2, &vars_vcf2))) {
         LOG_FATAL("%s\n", "Couldn't parse vars from 2nd vcf file");
         return -1;
    }
#if FIXME_TESTING
    for (i=0 ; i<10; i++) {
         vcf_write_var(vcfset_conf.vcf_out, vars_vcf2[i]);
    }
#endif

   vcf_parse_vars() vs parse_var()
   
   

         if (vcf_var_has_info_key(NULL, vcfset_conf.var, "INDEL")) {
         } else if (vcf_var_filtered(vcfset_conf.var)) {
    }
*/


    /* cleanup
     */

    for (i=0; i<num_vars_vcf2; i++) {
         vcf_free_var(& vars_vcf2[i]);
    }
    free(vars_vcf2);

    fclose(vcfset_conf.vcf_in1);
    fclose(vcfset_conf.vcf_in2);
    if (stdout != vcfset_conf.vcf_out) {
         fclose(vcfset_conf.vcf_out);
    }

    if (0==rc) {
         LOG_VERBOSE("%s\n", "Successful exit.");
    }

    LOG_ERROR("Unfinished\n");
    LOG_FIXME("1. Test me against lofreq_vcfset.py\n");
    LOG_FIXME("2. run valgrind\n");

    return rc;
}
/* main_vcfset */

