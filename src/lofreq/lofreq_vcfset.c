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



#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <assert.h>
#include <stdarg.h>
#include <getopt.h>
#include <stdlib.h>

#include "htslib/kstring.h"
#include "htslib/tbx.h"

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
     SETOP_CONCAT
} vcfset_op_t;

typedef struct {
     vcf_file_t vcf_in1;
     vcf_file_t vcf_in2;
     vcf_file_t vcf_out;
     vcfset_op_t vcf_setop;
     int only_passed; /* if 1, ignore any filtered variant */
     int only_pos; /* 0: allele aware. if 1, ignore ref and alt base during comparisons.  */
     int only_snvs;
     int only_indels;
} vcfset_conf_t;



static void
usage(const vcfset_conf_t* vcfset_conf)
{
     fprintf(stderr, "%s: Perform set operations on two vcf files\n\n", MYNAME);
     fprintf(stderr, "Usage: %s [options] -a op -1 1.vcf -2 2.vcf \n", MYNAME);

     fprintf(stderr,"Options:\n");
     fprintf(stderr, "  -1 | --vcf1 FILE      1st VCF input file (bgzip supported)\n");
     fprintf(stderr, "  -2 | --vcf2 FILE      2nd VCF input file (mandatory - except for concat - and needs to be tabix indexed)\n");
     fprintf(stderr, "  -o | --vcfout         VCF output file (default: - for stdout; gzip supported).\n");
     fprintf(stderr, "  -a | --action         Set operation to perform: intersect, complement or concat.\n"
             "                        - intersect = vcf1 AND vcf2.\n"
             "                        - complement = vcf1 \\ vcf2.\n"
             "                        - concat = vcf1 + vcf2 ... vcfn (output as in file order, i.e. output not necessarily sorted!)\n");
     fprintf(stderr, "  -I | --add-info STR   Add info field, e.g. 'SOMATIC'\n");
     fprintf(stderr, "       --count-only     Don't print bases, just numbers\n");
     fprintf(stderr, "       --only-pos       Disable allele-awareness by using position only (ignoring bases) as key for storing and comparison\n");
     fprintf(stderr, "       --only-passed    Ignore variants marked as filtered\n");
     fprintf(stderr, "       --only-snvs      Ignore anything but SNVs in both input files\n");
     fprintf(stderr, "       --only-indels    Ignore anything but indels in both input files\n");
     fprintf(stderr, "       --verbose        Be verbose\n");
     fprintf(stderr, "       --debug          Enable debugging\n");

     fprintf(stderr, "\nNote, vcf1 is always fully parsed, whereas indexing is used for vcf2.\n");
     fprintf(stderr, "Therefore, use the bigger file as vcf2 to speed things up.\n");
     fprintf(stderr, "Header/meta-data for the output file is taken from vcf1\n");
}
/* usage() */




int 
main_vcfset(int argc, char *argv[])
{
     vcfset_conf_t vcfset_conf;
     char *vcf_header = NULL;
     int rc = 0;
     char *vcf_in1, *vcf_in2, *vcf_out;
     long int num_vars_vcf1;
     long int num_vars_vcf1_ign, num_vars_out;
     static int only_passed = 0;
     static int only_pos = 0;
     static int only_snvs = 0;
     static int only_indels = 0;
     static int count_only = 0;
     tbx_t *vcf2_tbx = NULL; /* index for second vcf file */
     htsFile *vcf2_hts = NULL;
     char *add_info_field = NULL;
     int vcf_concat_findex = 0;
     vcf_in1 = vcf_in2 = vcf_out = NULL;
     num_vars_vcf1 = 0;
     num_vars_vcf1_ign = num_vars_out = 0;

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
              {"only-pos", no_argument, &only_pos, 1},
              {"only-indels", no_argument, &only_indels, 1},
              {"only-snvs", no_argument, &only_snvs, 1},
              {"count-only", no_argument, &count_only, 1},

              {"vcf1", required_argument, NULL, '1'},
              {"vcf2", required_argument, NULL, '2'},
              {"vcfout", required_argument, NULL, 'o'},
              {"action", required_argument, NULL, 'a'},
              {"add-info", required_argument, NULL, 'I'},

              {0, 0, 0, 0} /* sentinel */
         };

         /* keep in sync with long_opts and usage */
         static const char *long_opts_str = "h1:2:o:a:I:";

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

              } else if (0 == strcmp(optarg, "concat")) {
                   vcfset_conf.vcf_setop = SETOP_CONCAT;

              } else {
                   LOG_FATAL("Unknown action '%s'. Exiting...\n", optarg);
                   free(vcf_in1); free(vcf_in2); free(vcf_out);
                   return 1;
              }
              break;

         case 'I': 
              add_info_field = strdup(optarg);
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
    vcfset_conf.only_pos = only_pos;
    vcfset_conf.only_snvs = only_snvs;
    vcfset_conf.only_indels = only_indels;

    if (vcfset_conf.only_indels && vcfset_conf.only_snvs) {
         LOG_FATAL("%s\n", "Can't take only indels *and* only snvs into account");
         return 1;
    }

    if (0 != argc - optind - 1) {
         if (vcfset_conf.vcf_setop == SETOP_CONCAT) {
              vcf_concat_findex = optind;
         } else {
              LOG_FATAL("%s\n", "Unrecognized arguments found\n");
              return 1;
         }
    } else {
         if (vcfset_conf.vcf_setop == SETOP_CONCAT) {
              LOG_FATAL("%s\n", "No extra files for concat given\n");
              return 1;
         }
    }
#if 0
    int i; for (i=optind+1; i<argc; i++) {
         LOG_FIXME("argv[%d]=%s\n", i, argv[i]);
    }
#endif

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

    if  (vcf_in1 == NULL || (vcf_in2 == NULL && vcfset_conf.vcf_setop != SETOP_CONCAT)) {
         LOG_FATAL("%s\n\n", "At least one vcf input file not specified");
         usage(& vcfset_conf);
         free(vcf_in1); free(vcf_in2); free(vcf_out);
         return 1;
    }
    if (vcf_in2 != NULL && vcfset_conf.vcf_setop == SETOP_CONCAT) {
         LOG_FATAL("%s\n\n", "For concat just use the -1 option followed by all other vcf files instead of using -2");
         usage(& vcfset_conf);
         free(vcf_in1); free(vcf_in2); free(vcf_out);
         return 1;         
    }

    if (vcf_file_open(& vcfset_conf.vcf_in1, vcf_in1, 
                      HAS_GZIP_EXT(vcf_in1), 'r')) {
         LOG_ERROR("Couldn't open %s\n", vcf_in1);
         free(vcf_in1); free(vcf_in2); free(vcf_out);
         return 1;
    }

    if (vcf_in2) {
         vcf2_hts = hts_open(vcf_in2, "r");
         if (!vcf2_hts) {
              LOG_FATAL("Couldn't load %s\n", vcf_in2);
              return 1;
         }
         vcf2_tbx = tbx_index_load(vcf_in2);
         if (!vcf2_tbx) {
              LOG_FATAL("Couldn't load tabix index for %s\n", vcf_in2);
              return 1;
         }
    }

    /* vcf_out default if not set: stdout==- */
    if (! vcf_out) {
         vcf_out = malloc(2 * sizeof(char));
         strcpy(vcf_out, "-");
    }

    if (! count_only) {
         if (vcf_file_open(& vcfset_conf.vcf_out, vcf_out, 
                           HAS_GZIP_EXT(vcf_out), 'w')) {
              LOG_ERROR("Couldn't open %s\n", vcf_out);
              free(vcf_in1); free(vcf_in2); free(vcf_out);
              return 1;
         }
    }

    /* use meta-data/header of vcf_in1 for output
     */
    LOG_DEBUG("Getting header from %s\n", vcf_in1);
    if (0 !=  vcf_parse_header(&vcf_header, & vcfset_conf.vcf_in1)) {
         LOG_WARN("%s\n", "vcf_parse_header() failed");
         if (vcf_file_seek(& vcfset_conf.vcf_in1, 0, SEEK_SET)) {
              LOG_FATAL("%s\n", "Couldn't rewind file to parse variants"
                        " after header parsing failed");
              return -1;
         }
    } else {
         if (! count_only) {
              /* vcf_write_header would write *default* header */
              vcf_write_header(& vcfset_conf.vcf_out, vcf_header);
         }
         free(vcf_header);
    }

    
    /* parse first vcf file
     */
    LOG_DEBUG("Starting to parse variants from %s\n", vcf_in1);
    while (1) {
         var_t *var1 = NULL;
         int rc;
         int is_indel;
         kstring_t var2_kstr = {0, 0, 0};
         hts_itr_t *var2_itr = NULL;
         char regbuf[1024];
         int var2_match = 0;

         vcf_new_var(&var1);
         rc = vcf_parse_var(& vcfset_conf.vcf_in1, var1);
         if (rc) {
              free(var1);
              
              if (vcfset_conf.vcf_setop != SETOP_CONCAT) {
                   break;
              } else {
                   vcf_concat_findex++;
                   if (vcf_concat_findex==argc) {
                        break;
                   }
                   /* set vcf1 up anew and simply continue as if nothing happened 
                    */
                   vcf_file_close(& vcfset_conf.vcf_in1);
                   free(vcf_in1);

                   vcf_in1 = strdup(argv[vcf_concat_findex]);
                   LOG_DEBUG("updated vcf_in1 = %s\n", vcf_in1);
                   if (vcf_file_open(& vcfset_conf.vcf_in1, vcf_in1, 
                                     HAS_GZIP_EXT(vcf_in1), 'r')) {
                        LOG_ERROR("Couldn't open %s\n", vcf_in1);
                        free(vcf_in1); free(vcf_in2); free(vcf_out);
                        return 1;
                   }
                   if (0 != vcf_skip_header(& vcfset_conf.vcf_in1)) {
                        LOG_WARN("skip header failed for %s\n", vcf_in1);
                   }
                   continue;
              }
         }

         is_indel = vcf_var_is_indel(var1);
         if (vcfset_conf.only_snvs && is_indel) {
              free(var1);
              continue;
         } else if (vcfset_conf.only_indels && ! is_indel) {
              free(var1);
              continue;
         }

         if (! vcfset_conf.only_pos && NULL != strchr(var1->alt, ',')) {
              LOG_FATAL("%s\n", "No support for multi-allelic SNVs in vcf1");
              return -1;
         }
         if (vcfset_conf.only_passed && ! VCF_VAR_PASSES(var1)) {
#ifdef TRACE
              LOG_DEBUG("Skipping non-passing var1 %s:%d\n", var1->chrom, var1->pos);
#endif
              num_vars_vcf1_ign += 1;
              vcf_free_var(& var1);
              continue;
         }
         if (add_info_field) {
              vcf_var_add_to_info(var1, add_info_field);
         }
         num_vars_vcf1 += 1;
#ifdef TRACE
         LOG_DEBUG("Got passing var1 %s:%d\n", var1->chrom, var1->pos);
#endif

         if (vcfset_conf.vcf_setop == SETOP_CONCAT) {
              num_vars_out += 1;
              if (! count_only) {
                   vcf_write_var(& vcfset_conf.vcf_out, var1);
              }
              vcf_free_var(& var1);
              /* skip comparison against vcf2 */
              continue;
         }

         /* use index access to vcf2 */
         snprintf(regbuf, 1024, "%s:%ld-%ld", var1->chrom, var1->pos+1, var1->pos+1);
         var2_itr = tbx_itr_querys(vcf2_tbx, regbuf);
         if (! var2_itr) {
              var2_match = 0;
         } else {
              var2_match = 0;
              while (tbx_itr_next(vcf2_hts, vcf2_tbx, var2_itr, &var2_kstr) >= 0) {
                   var_t *var2 = NULL;
                   int var2_is_indel = 0;

                   vcf_new_var(&var2);
                   rc = vcf_parse_var_from_line(var2_kstr.s, var2);
                   /* LOG_FIXME("%d:%s>%s looking at var2 %d:%s>%s (reg %s)\n", 
                             var1->pos+1, var1->ref, var1->alt,
                             var2->pos+1, var2->ref, var2->alt, regbuf); */
                   if (rc) {
                        LOG_FATAL("%s\n", "Error while parsing variant returned from tabix");
                        return -1;
                   }

                   var2_is_indel = vcf_var_is_indel(var2);

                   /* iterator returns anything overlapping with that 
                    * position, i.e. this also includes up/downstream
                    * indels, so make sure actual position matches */
                   if (var1->pos != var2->pos) {
                        var2_match = 0;

                   } else if (vcfset_conf.only_passed && ! VCF_VAR_PASSES(var2)) {
                        var2_match = 0;

                   } else if (vcfset_conf.only_snvs && var2_is_indel) {
                        var2_match = 0;

                   } else if (vcfset_conf.only_indels && ! var2_is_indel) {
                        var2_match = 0;

                   } else if (vcfset_conf.only_pos) {
#ifdef TRACE
                        LOG_DEBUG("Pos match for var2 %s:%d\n", var2->chrom, var2->pos);
#endif
                        var2_match = 1;

                   } else {
                        if (0==strcmp(var1->ref, var2->ref) && 0==strcmp(var1->alt, var2->alt)) {
#ifdef TRACE
                             LOG_DEBUG("Full match for var2 %s:%d\n", var2->chrom, var2->pos);
#endif
                             var2_match = 1;/* FIXME: check type as well i.e. snv vs indel */                             
                        }
                   }
                   vcf_free_var(&var2);
                   if (var2_match) {
                        break;/* no need to continue */
                   }
              }
         }

         if (vcfset_conf.vcf_setop == SETOP_COMPLEMENT) {
              /* relative complement : elements in A but not B */
              if (!var2_match) {
                   num_vars_out += 1;
                   if (! count_only) {
                        vcf_write_var(& vcfset_conf.vcf_out, var1);
                   }
              }
         } else if (vcfset_conf.vcf_setop == SETOP_INTERSECT) {
              if (var2_match) {
                   num_vars_out += 1;
                   if (! count_only) {
                        vcf_write_var(& vcfset_conf.vcf_out, var1);
                   }
              }

         } else {
              LOG_FATAL("Internal error: unsupported vcf_setop %d\n", vcfset_conf.vcf_setop);
              return 1;
         }

         vcf_free_var(& var1);
         tbx_itr_destroy(var2_itr);
    }/* while (1) */

    vcf_file_close(& vcfset_conf.vcf_in1);
    if (vcf_in2) {
         hts_close(vcf2_hts);
         tbx_destroy(vcf2_tbx);
    }
    LOG_VERBOSE("Parsed %d variants from 1st vcf file (ignoring %d non-passed of those)\n", 
                num_vars_vcf1 + num_vars_vcf1_ign, num_vars_vcf1_ign);
    LOG_VERBOSE("Wrote %d variants to output\n", 
                num_vars_out);
    if (! count_only) {
         vcf_file_close(& vcfset_conf.vcf_out);
    }

    if (0==rc) {
         if (count_only) {
              printf("%ld\n", num_vars_out);
         }

         LOG_VERBOSE("%s\n", "Successful exit.");
    }

    free(vcf_in1);
    free(vcf_in2);
    free(vcf_out);


    return rc;
}
/* main_vcfset */

