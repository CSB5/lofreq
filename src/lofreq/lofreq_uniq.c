/* -*- c-file-style: "k&r"; indent-tabs-mode: nil; -*- */


/*********************************************************************
 *
 * FIXME update license
 *
 *********************************************************************/


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

/* lofreq includes */
#include "vcf.h"
#include "utils.h"
#include "log.h"
#include "binom.h"
#include "plp.h"



#define DEFAULT_SIG 0.05

#define BUF_SIZE 1<<16

typedef struct {
     /* fix values once established */
     float sig;
     FILE *vcf_out;
     FILE *vcf_in;
     long long int bonf;

     /* changing per pos: the var to test */
     var_t *var;
} uniq_conf_t;



/* used as pileup callback function which is not ideal since this can
 * only work on one position (has to be ensured by caller).
 *
 * No cov means I won't be called and no output will be generated.
 * Non-sig pv means I'm not sure and no ouput will be generated. Only
 * if pv is sig we will print the var
 */
void
uniq_snv(const plp_col_t *p, void *confp)
{
     uniq_conf_t *conf = (uniq_conf_t *)confp;
     int alt_count;
     char *af_char = NULL;
     float af;
     double pvalue = DBL_MAX;
     int is_sig;

     if (vcf_var_has_info_key(NULL, conf->var, "INDEL")) {
          LOG_WARN("uniq logic can't be applied to indels."
                   " Skipping indel var at %s %d\n", 
                   conf->var->chrom, conf->var->pos+1);
          return;
     }

     if (0 != strcmp(p->target, conf->var->chrom) || p->pos != conf->var->pos) {
          LOG_ERROR("wrong pileup for var. pileup for %s %d. var for %s %d\n", 
                    p->target, p->pos+1, conf->var->chrom, conf->var->pos+1);
          return;
     }

     /* no cov usually means I won't be called, but just to be safe: 
      */
     if (0 == p->coverage) {          
          return;
     }

     vcf_var_has_info_key(&af_char, conf->var, "AF");
     if (NULL == af_char) {
          LOG_ERROR("%s\n", "Couldn't parse AF (key not found) from the following variant:");
          vcf_write_var(stderr, conf->var);
          return;
     }
     af = strtof(af_char, (char **)NULL); /* atof */
     free(af_char);
     if (af < 0.0 || af > 1.0) {
          LOG_ERROR("%s\n", "Couldn't parse AF (value out of bound) from the following variant:");
          vcf_write_var(stderr, conf->var);
          return;
     }

     alt_count = base_count(p, conf->var->alt);

#ifdef DEBUG
     LOG_DEBUG("Now testing af=%f cov=%d alt_count=%d at %s %d for var:",
               af, p->coverage, alt_count, p->target, p->pos+1);
#endif

     if (0 != binom_sf(&pvalue, p->coverage, alt_count, af)) {
          LOG_ERROR("%s\n", "binom_sf() failed");
          return;
     }

     is_sig = pvalue < conf->sig/(float)conf->bonf;
     LOG_VERBOSE("%s %d %c>%c AF=%f | %s (p-value=%g sig/bonf=%g) | BAM alt_count=%d cov=%d (freq=%f)\n", 
                 conf->var->chrom, conf->var->pos+1, conf->var->ref, conf->var->alt, af,
                 is_sig ? "unique" : "not necessarily unique", pvalue, conf->sig/(float)conf->bonf,
                 alt_count, p->coverage, alt_count/(float)p->coverage);
     if (is_sig) {
          vcf_write_var(stderr, conf->var);
     }
}


static void
usage(const uniq_conf_t* uniq_conf)
{
     fprintf(stderr, "Usage: %s uniq [options] in.bam\n\n", PACKAGE);
     fprintf(stderr, "Options:\n");
     /* generic */
     fprintf(stderr, "       --verbose            Be verbose\n");
     fprintf(stderr, "       --debug              Enable debugging\n");
     fprintf(stderr, "  -v | --vcf-in FILE       Input vcf file listing variants [- = stdin]\n");
     fprintf(stderr, "  -o | --vcf-out FILE      Output vcf file [- = stdout]\n");
     fprintf(stderr, "  -s | --sig               Significance threshold [%f]\n", uniq_conf->sig);
}
/* usage() */



int 
main_uniq(int argc, char *argv[])
{
     int c, i;
     char *bam_file = NULL;
     char *vcf_in = NULL; /* - == stdout */
     char *vcf_out = NULL; /* - == stdout */
     mplp_conf_t mplp_conf;
     uniq_conf_t uniq_conf;
     void (*plp_proc_func)(const plp_col_t*, void*);
     int rc = 0;
     var_t **vars = NULL;
     int num_vars = 0;
     char *vcf_header;

     for (i=0; i<argc; i++) {
          LOG_DEBUG("arg %d: %s\n", i, argv[i]);
     }


     /* default uniq options */
     memset(&uniq_conf, 0, sizeof(uniq_conf_t));
     uniq_conf.vcf_in = stdin;
     uniq_conf.vcf_out = stdout;
     uniq_conf.sig = DEFAULT_SIG;
     uniq_conf.bonf = 1;

     /* default pileup options */
     memset(&mplp_conf, 0, sizeof(mplp_conf_t));
     mplp_conf.max_mq = 255;
     mplp_conf.min_bq = 3;
     mplp_conf.capQ_thres = 0;
     mplp_conf.max_depth = 1000000;
     mplp_conf.flag = MPLP_NO_ORPHAN;
    

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

              {"vcf-in", required_argument, NULL, 'v'},
              {"vcf-out", required_argument, NULL, 'o'},
              {"sig", required_argument, NULL, 's'},

              {0, 0, 0, 0} /* sentinel */
         };

         /* keep in sync with long_opts and usage */
         static const char *long_opts_str = "hv:o:s:"; 

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
              usage(& uniq_conf); 
              return 0;

         case 'v': 
              if (0 != strcmp(optarg, "-")) {
                   if (! file_exists(optarg)) {
                        LOG_FATAL("Input file '%s' does not exist. Exiting...\n", optarg);
                        return 1;
                   }
              }
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

         case 's': 
              uniq_conf.sig = strtof(optarg, (char **)NULL); /* atof */
              if (0==uniq_conf.sig) {
                   LOG_FATAL("%s\n", "Couldn't parse sign-threshold"); 
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
        usage(& uniq_conf);
        return 1;
    }

    if (1 != argc - optind - 1) {
        fprintf(stderr, "Need exactly one BAM file as last argument\n");
        return 1;
    }
    bam_file = (argv + optind + 1)[0];
    if (! file_exists(bam_file)) {
         LOG_FATAL("BAM file %s does not exist. Exiting...\n", bam_file);
         return -1;
    }

    if (NULL != vcf_out && 0 != strcmp(vcf_out, "-")) {
         uniq_conf.vcf_out = fopen(vcf_out, "w");
    }

    if (NULL != vcf_in && 0 != strcmp(vcf_in, "-")) {
         uniq_conf.vcf_in = fopen(vcf_in, "r");
    }
    if (0 !=  vcf_parse_header(&vcf_header, uniq_conf.vcf_in)) {
         LOG_FATAL("%s\n", "vcf_parse_header() failed");
         free(vcf_header);
         return -1;
    } 
    fprintf(uniq_conf.vcf_out, "%s", vcf_header);
    free(vcf_header);

    if (-1 == (num_vars = vcf_parse_vars(uniq_conf.vcf_in, &vars))) {
         LOG_FATAL("%s\n", "vcf_parse_vars() failed");
         return -1;
    }

    plp_proc_func = &uniq_snv;
    uniq_conf.bonf = num_vars;
    for (i=0; i<num_vars; i++) {
         char reg_buf[BUF_SIZE];

         uniq_conf.var = vars[i];

         snprintf(reg_buf, BUF_SIZE, "%s:%ld-%ld", 
                  vars[i]->chrom, vars[i]->pos+1, vars[i]->pos+1);
         mplp_conf.reg = strdup(reg_buf);

         LOG_DEBUG("pileup for var no %d at %s %d\n", 
                   i+1, uniq_conf.var->chrom, uniq_conf.var->pos+1);
         if (vcf_var_has_info_key(NULL, uniq_conf.var, "INDEL")) {
              LOG_WARN("Skipping indel var at %s %d\n", 
                       uniq_conf.var->chrom, uniq_conf.var->pos+1);
              continue;
         }

         rc = mpileup(&mplp_conf, plp_proc_func, (void*)&uniq_conf,
                      1, (const char **) argv + optind + 1);

         free(mplp_conf.reg); 
         mplp_conf.reg = NULL;
    }

    if (stdin != uniq_conf.vcf_in) {
         fclose(uniq_conf.vcf_in);
    }
    if (stdout != uniq_conf.vcf_out) {
         fclose(uniq_conf.vcf_out);
    }

    for (i=0; i<num_vars; i++) {
         vcf_free_var(& vars[i]);
    }
    free(vars);

    free(vcf_in);
    free(vcf_out);

    if (0==rc) {
         LOG_VERBOSE("%s\n", "Successful exit.");
    }

    LOG_FIXME("%s\n", "test: indel skip");
    LOG_FIXME("%s\n", "add: ign-filtered/passed-only");
    LOG_FIXME("%s\n", "test: against self bam: denv2-simulation/denv2-10haplo_true-snp.vcf denv2-simulation/denv2-10haplo.bam should produce no SNVs");
    LOG_FIXME("%s\n", "test: region access fast?");

    return rc;
}
/* main_call */

