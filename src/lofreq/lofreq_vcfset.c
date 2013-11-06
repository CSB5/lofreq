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
     vcf_file_t vcf_in1;
     vcf_file_t vcf_in2;
     vcf_file_t vcf_out;
     vcfset_op_t vcf_setop;
     int only_passed; /* if 1, ignore any filtered variant */
     int use_bases; /* if 0, ignore ref and alt base during comparisons */
} vcfset_conf_t;



static void
usage(const vcfset_conf_t* vcfset_conf)
{
     fprintf(stderr, "%s: Perform set operations on two vcf files\n\n", MYNAME);
     fprintf(stderr, "Usage: %s [options] -a op -1 1.vcf -2 2.vcf \n", MYNAME);

     fprintf(stderr,"Options:\n");
     fprintf(stderr, "  -1 | --vcf1 FILE      1st VCF input file (gzip supported)\n");
     fprintf(stderr, "  -2 | --vcf2 FILE      2nd VCF input file (gzip supported)\n");
     fprintf(stderr, "  -o | --vcfout         VCF output file (default: - for stdout; gzip supported).\n"
             "                        Meta-data will be copied from vcf1\n");
     fprintf(stderr, "  -a | --action         Set operation to perform:\n"
             "                        intersect or complement.\n"
             "                        intersect = vcf1 AND vcf2.\n"
/*             "                        union = vcf1 OR vcf2.\n"*/
             "                        complement = vcf1 \\ vcf2.\n");
     fprintf(stderr, "       --use-bases      Use position and bases as key for storing and comparison. By default only position is used\n");
     fprintf(stderr, "       --only-passed    Ignore variants marked as filtered\n");
     fprintf(stderr, "       --verbose        Be verbose\n");
     fprintf(stderr, "       --debug          Enable debugging\n");
}
/* usage() */




int 
main_vcfset(int argc, char *argv[])
{
     vcfset_conf_t vcfset_conf;
     char *vcf_header = NULL;
     int rc = 0;
     char *vcf_in1, *vcf_in2, *vcf_out;
     long int num_vars_vcf1, num_vars_vcf2;
     long int num_vars_vcf1_ign, num_vars_vcf2_ign, num_vars_out;
     var_hash_t *var_hash_vcf2 = NULL; /* must be declared NULL ! */
     static int only_passed = 0;
     static int use_bases = 0;

     vcf_in1 = vcf_in2 = vcf_out = NULL;
     num_vars_vcf1 = num_vars_vcf2 = 0;
     num_vars_vcf1_ign = num_vars_vcf2_ign = num_vars_out = 0;

     /* default vcfset options */
     memset(&vcfset_conf, 0, sizeof(vcfset_conf_t));
     /* vcfset_conf.vcf_in1 = NULL; */
     /* vcfset_conf.vcf_in2 = NULL; */
     /* vcfset_conf.vcf_out = stdout;*/


    /* keep in sync with long_opts_str and usage 
     *
     * getopt is a pain in the whole when it comes to syncing of long
     * and short args and usage. check out gopt, libcfu...
     */
    while (1) {
         int c;
         static struct option long_opts[] = {
              /* see usage sync */
              {"help", no_argument, NULL, 'h'},
              {"verbose", no_argument, &verbose, 1},
              {"debug", no_argument, &debug, 1},
              {"only-passed", no_argument, &only_passed, 1},
              {"use-bases", no_argument, &use_bases, 1},

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
              free(vcf_in1); free(vcf_in2); free(vcf_out);
              return 0;

         case '1': 
              vcf_in1 = strdup(optarg);
              break;

         case '2': 
              vcf_in2 = strdup(optarg);
              break;

         case 'o':
              if (0 != strcmp(optarg, "-")) {
                   if (file_exists(optarg)) {
                        LOG_FATAL("Cowardly refusing to overwrite file '%s'. Exiting...\n", optarg);
                        free(vcf_in1); free(vcf_in2);
                        return 1;
                   }
              }
              vcf_out = strdup(optarg);
              break;

         case 'a': 
              if (0 == strcmp(optarg, "intersect")) {
                   vcfset_conf.vcf_setop = SETOP_INTERSECT;

              } else if (0 == strcmp(optarg, "complement")) {
                   vcfset_conf.vcf_setop = SETOP_COMPLEMENT;

/*              } else if (0 == strcmp(optarg, "union")) {
                vcfset_conf.vcf_setop = SETOP_UNION; */

              } else {
                   LOG_FATAL("Unknown action '%s'. Exiting...\n", optarg);
                   free(vcf_in1); free(vcf_in2); free(vcf_out);
                   return 1;
              }
              break;

         case '?': 
              LOG_FATAL("%s\n", "unrecognized arguments found. Exiting...\n"); 
              free(vcf_in1); free(vcf_in2); free(vcf_out);
              return 1;

         default:
              break;
         }
    }

    vcfset_conf.only_passed = only_passed;
    vcfset_conf.use_bases = use_bases;


    if (0 != argc - optind - 1) {
         LOG_FATAL("%s\n", "Unrecognized arguments found\n");
         return 1;
    }

    if (argc == 2) {
        fprintf(stderr, "\n");
        usage(& vcfset_conf);
        free(vcf_in1); free(vcf_in2); free(vcf_out);
        return 1;
    }

    if (vcfset_conf.vcf_setop == SETOP_UNKNOWN) {
         LOG_FATAL("%s\n", "No set operation specified");
         usage(& vcfset_conf);
         free(vcf_in1); free(vcf_in2); free(vcf_out);
         return 1;
    }

    /* all input files must be specified */
    if  (vcf_in1 == NULL || vcf_in2 == NULL) {
         LOG_FATAL("%s\n\n", "At least one needed the vcf input file not specified");
         usage(& vcfset_conf);
         free(vcf_in1); free(vcf_in2); free(vcf_out);
         return 1;
         
    }

    /* set up all input files */
    if (vcf_file_open(& vcfset_conf.vcf_in1, vcf_in1, 
                      HAS_GZIP_EXT(vcf_in1), 'r')) {
         LOG_ERROR("Couldn't open %s\n", vcf_in1);
         free(vcf_in1); free(vcf_in2); free(vcf_out);
         return 1;
    }
    if (vcf_file_open(& vcfset_conf.vcf_in2, vcf_in2, 
                      HAS_GZIP_EXT(vcf_in2), 'r')) {
         LOG_ERROR("Couldn't open %s\n", vcf_in2);
         free(vcf_in1); free(vcf_in2); free(vcf_out);
         return 1;
    }


    /* vcf_out default if not set: stdout==- */
    if (! vcf_out) {
         vcf_out = malloc(2 * sizeof(char));
         strcpy(vcf_out, "-");
    }

    if (vcf_file_open(& vcfset_conf.vcf_out, vcf_out, 
                      HAS_GZIP_EXT(vcf_out), 'w')) {
         LOG_ERROR("Couldn't open %s\n", vcf_out);
         free(vcf_in1); free(vcf_in2); free(vcf_out);
         return 1;
    }

    free(vcf_in1);
    free(vcf_in2);
    free(vcf_out);


    /* recipe: read B into memory and parse from A one by one
     * ======================================================
     *
     * things can be done a lot more efficient if both input files are
     * sorted, but let's assume they are not.
     *
     */

    /* use meta-data/header of vcf_in1 for output
     */
    if (0 !=  vcf_parse_header(&vcf_header, & vcfset_conf.vcf_in1)) {
         LOG_WARN("%s\n", "vcf_parse_header() failed");
         if (vcf_file_seek(& vcfset_conf.vcf_in1, 0, SEEK_SET)) {
              LOG_FATAL("%s\n", "Couldn't rewind file to parse variants"
                        " after header parsing failed");
              return -1;
         }
    } else {
         /* vcf_write_header would write *default* header */
         vcf_write_header(& vcfset_conf.vcf_out, vcf_header);
         free(vcf_header);
    }


    /* skip meta-data/header in vcf_in2
     */
    if (0 != vcf_skip_header(& vcfset_conf.vcf_in2)) {
         LOG_FATAL("%s\n", "Failed to skip header in 2nd vcf files");
         return -1;
    }


    /* read vcf2 and save as hash
     */
    while (1) {
         var_t *var;
         char *key;
         int rc;
         vcf_new_var(&var);
         rc = vcf_parse_var(& vcfset_conf.vcf_in2, var);
         if (-1 == rc) {
              LOG_FATAL("%s\n", "Parsing error while parsing 2nd vcf-file");
              exit(1);
         }
         if (1 == rc) {/* EOF */
              free(var);
              break;
         }

         if (! vcfset_conf.only_passed || VCF_VAR_PASSES(var)) {
              vcfset_conf.use_bases ? vcf_var_key(&key, var) : vcf_var_key_simple(&key, var);
              var_hash_add(& var_hash_vcf2, key, var);
#ifdef TRACE
              fprintf(stderr, "var_2 pass: "); vcf_write_var(stderr, var);
#endif
         } else {
              num_vars_vcf2_ign += 1;
              vcf_free_var(& var);
         }
#ifdef TRACE
         LOG_DEBUG("Adding %s\n", key);
#endif
    }
    vcf_file_close(& vcfset_conf.vcf_in2);

    num_vars_vcf2 = HASH_COUNT(var_hash_vcf2);
    LOG_VERBOSE("Parsed %d variants from 2nd vcf file (ignoring %d non-passed of those)\n", 
                num_vars_vcf2 + num_vars_vcf2_ign, num_vars_vcf2_ign);


    /* now parse first vcf file and decide what to do
     */
    while (1) {
         var_t *var_1;
         char *key;
         var_hash_t *var_2;
         int rc;

         vcf_new_var(&var_1);
         rc = vcf_parse_var(& vcfset_conf.vcf_in1, var_1);
         if (-1 == rc) {
              LOG_FATAL("%s\n", "Parsing error while parsing 1st vcf-file");
              exit(1);
         }
         if (1 == rc) {/* EOF */
              free(var_1);
              break;
         }

         if (vcfset_conf.only_passed && ! VCF_VAR_PASSES(var_1)) {
              num_vars_vcf1_ign += 1;
              vcf_free_var(& var_1);
              continue;
         }
#ifdef TRACE
         fprintf(stderr, "var_1 pass: "); vcf_write_var(stderr, var_1);
#endif
         num_vars_vcf1 += 1;
         vcfset_conf.use_bases ? vcf_var_key(&key, var_1) : vcf_var_key_simple(&key, var_1);
         HASH_FIND_STR(var_hash_vcf2, key, var_2);
#ifdef TRACE
         LOG_DEBUG("var with key %s in 2: %s\n", key, var_2? "found" : "not found");
#endif
         free(key);

         if (vcfset_conf.vcf_setop == SETOP_COMPLEMENT) {
              /* relative complement : elements in A but not B */
              if (NULL == var_2) {
                   num_vars_out += 1;
                   vcf_write_var(& vcfset_conf.vcf_out, var_1);

              }
#ifdef NOT_OKAY_IF_SIMPLE_KEY_AND_MULTIALLELIC
              else {
                   /* save some mem */
                   HASH_DEL(var_hash_vcf2, var_2);
                   var_hash_free_elem(var_2);
              }
#endif
         } else if (vcfset_conf.vcf_setop == SETOP_INTERSECT) {
              if (NULL != var_2) {
                   num_vars_out += 1;
                   vcf_write_var(& vcfset_conf.vcf_out, var_1);
              }

         } else {
              LOG_FATAL("Internal error: unsupported vcf_setop %d\n", vcfset_conf.vcf_setop);
              return 1;
         }

         vcf_free_var(& var_1);
    }
    vcf_file_close(& vcfset_conf.vcf_in1);


    LOG_VERBOSE("Parsed %d variants from 1st vcf file (ignoring %d non-passed of those)\n", 
                num_vars_vcf1 + num_vars_vcf1_ign, num_vars_vcf1_ign);
    LOG_VERBOSE("Wrote %d variants to output\n", 
                num_vars_out);
    vcf_file_close(& vcfset_conf.vcf_out);

    var_hash_free_table(var_hash_vcf2);

    if (0==rc) {
         LOG_VERBOSE("%s\n", "Successful exit.");
    }

    return rc;
}
/* main_vcfset */

