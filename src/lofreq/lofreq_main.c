/* -*- c-file-style: "k&r"; indent-tabs-mode: nil; -*- */

/*********************************************************************
 *
 * FIXME update license
 *
 *********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <libgen.h>

/* lofreq includes */
#include "log.h"
#include "utils.h"
#ifdef USE_ALNERRPROF
#include "lofreq_bamstats.h"
#endif
#include "lofreq_checkref.h"
#include "lofreq_filter.h"
#include "lofreq_snpcaller.h"
#include "lofreq_uniq.h"
#include "lofreq_vcfset.h"
#include "lofreq_index.h"

/*added lofreq_viterbi.h */
#include "lofreq_viterbi.h"

#ifndef __DATE__
__DATE__ = "NA";
#endif


static void prepend_dir_to_path(const char *dir_to_add)
{
     char *old_path = NULL;
     char *new_path;
     const char *PATH_NAME = "PATH";

     old_path = getenv(PATH_NAME);
     if (NULL == old_path) {
          setenv(PATH_NAME, dir_to_add, 1);
          return;
     }

     new_path = malloc((strlen(dir_to_add) + 1 + strlen(old_path) + 1) * sizeof(char));
     sprintf(new_path, "%s:%s", dir_to_add, old_path);
#if 0
     LOG_WARN("old path: %s\n", old_path);
     LOG_WARN("new path: %s\n", new_path);
#endif
     setenv(PATH_NAME, new_path, 1);
     free(new_path);
}



/* prepend dirname(argv0) and python source dir to PATH. This way we
 * make sure that package works even without properly installing it
 * and that the binary can repeatedly can call itself it necessary
 */
void 
add_local_dir_to_path(char *argv0) {
     char argv0_resolved[PATH_MAX];
     char *dirname_argv0 = NULL;
     char *lofreq_script_abs = NULL;
     char *lofreq_tools_abs = NULL;
     char lofreq_script_rel[] = "../scripts/lofreq2_somatic.py";
     char lofreq_tools_rel[] ="../tools/scripts/lofreq2_vcfplot.py";
     char *add_dir = NULL;

     /* add lofreq dir 
      */
     if (NULL == realpath(argv0, argv0_resolved)) {
          return;          
     }
     if (NULL == (dirname_argv0 = strdup(dirname(argv0_resolved)))) {
          return;
     }     
     prepend_dir_to_path(dirname_argv0);


     /* add local script dir if present
      */
     if (NULL == (lofreq_script_abs = strdup(dirname_argv0))) {
          free(dirname_argv0);
          return;
     }
     if (NULL == join_paths(&lofreq_script_abs, lofreq_script_rel)) {
#if 0
          LOG_WARN("join_paths %s and %s failed\n", lofreq_script_abs, lofreq_script_rel);
#endif
          free(lofreq_script_abs);
          free(dirname_argv0);
          return;          
     }
     if (! file_exists(lofreq_script_abs)) {
#if 0
          LOG_WARN("%s doesnt' exist\n", lofreq_script_abs);
#endif
          free(lofreq_script_abs);
          free(dirname_argv0);
          return;
     }
     add_dir = strdup(dirname(lofreq_script_abs));
     prepend_dir_to_path(add_dir);
     free(lofreq_script_abs);
     free(add_dir);


     /* add local tools dir if present
      */
     if (NULL == (lofreq_tools_abs = strdup(dirname_argv0))) {
          free(dirname_argv0);
          return;
     }
     if (NULL == join_paths(&lofreq_tools_abs, lofreq_tools_rel)) {
#if 0
          LOG_WARN("join_paths %s and %s failed\n", lofreq_tools_abs, lofreq_tools_rel);
#endif
          free(lofreq_tools_abs);
          free(dirname_argv0);
          return;          
     }
     if (! file_exists(lofreq_tools_abs)) {
#if 0
          LOG_WARN("%s doesnt' exist\n", lofreq_tools_abs);
#endif
          free(lofreq_tools_abs);
          free(dirname_argv0);
          return;
     }
     add_dir = strdup(dirname(lofreq_tools_abs));
     prepend_dir_to_path(add_dir);
     free(lofreq_tools_abs);
     free(add_dir);

     free(dirname_argv0);
}



static void usage(const char *myname)
{
     fprintf(stderr, "\n%s: Fast and sensitive inference of single-nucleotide variants\n", PACKAGE_NAME);
     /*fprintf(stderr, "Version %s\n", PACKAGE_VERSION);*/
     fprintf(stderr, "\n");
     fprintf(stderr, "Usage: %s <command> [options]\n\n", myname);
     fprintf(stderr, "  Main Commands:\n");
     fprintf(stderr, "    call          : Call variants\n");
#ifdef PARALLEL
     fprintf(stderr, "    call-parallel : Call variants in parallel\n");
#endif
     fprintf(stderr, "    somatic       : Call somatic variants\n");
     fprintf(stderr, "\n");
     fprintf(stderr, "  Other Commands:\n");
     fprintf(stderr, "    checkref      : Check that reference fasta and BAM file match\n");
     fprintf(stderr, "    filter        : Filter variants in VCF file\n");
     fprintf(stderr, "    uniq          : Test whether variants predicted in only one sample really are unique\n");
     fprintf(stderr, "    plpsummary    : Print pileup summary per position\n"); 
#ifdef USE_ALNERRPROF
     fprintf(stderr, "    bamstats      : Collect BAM statistics\n");
#endif
     fprintf(stderr, "    vcfset        : VCF set operations\n");
	/*added the viterbi printout */
     fprintf(stderr, "    viterbi       : viterbi realignment\n");

     fprintf(stderr, "    version       : Print version info\n");
     fprintf(stderr, "\n");
     fprintf(stderr, "  Extra Tools (if installed):\n");
     fprintf(stderr, "    vcfplot       : Plot VCF statistics\n");
     fprintf(stderr, "    cluster       : Cluster variants in VCF file (supports legacy SNP format)\n");
     fprintf(stderr, "\n");
     fprintf(stderr, "  Samtools Clones:\n");
     fprintf(stderr, "    index         : Create index for BAM file\n");
     fprintf(stderr, "    idxstats      : Print stats for indexed BAM file\n");
     fprintf(stderr, "\n");

     fprintf(stderr, "\n");
}




int main(int argc, char *argv[])
{
     add_local_dir_to_path(argv[0]);

     if (argc < 2) {
          usage(BASENAME(argv[0]));
          return 1;
     }
     if (strcmp(argv[1], "call") == 0)  {
          return main_call(argc, argv);

     } else if (strcmp(argv[1], "uniq") == 0)  {
          return main_uniq(argc, argv);

     } else if (strcmp(argv[1], "vcfset") == 0)  {
          return main_vcfset(argc, argv);
     } else if (strcmp(argv[1], "viterbi") == 0){
	  return main_viterbi(argc,argv);  /*added the viterbi call!!! */
     } else if (strcmp(argv[1], "index") == 0)  {
          return main_index(argc, argv);

     } else if (strcmp(argv[1], "idxstats") == 0)  {
          return main_idxstats(argc, argv);

     } else if (strcmp(argv[1], "checkref") == 0) {
          return main_checkref(argc, argv);
          
     } else if (strcmp(argv[1], "info") == 0) {
          LOG_FIXME("%s\n", "NOT IMPLEMENTED YET: has BI, has BD, readlen, has extra BAQ. is_paired. all based on first, say, 10k read\n");
          return 1;

     } else if (strcmp(argv[1], "wizard") == 0) {
          LOG_FIXME("%s\n", "NOT IMPLEMENTED YET\n");
          return 1;

     } else if (strcmp(argv[1], "filter") == 0) {
          return main_filter(argc, argv);

     } else if (strcmp(argv[1], "somatic") == 0 ||
                strcmp(argv[1], "vcfplot") == 0 ||
                strcmp(argv[1], "call-parallel") == 0 ||
                strcmp(argv[1], "cluster") == 0) {
          char **argv_execvp = calloc(argc, sizeof(char*));
          int i;
          char *somatic_script = "lofreq2_somatic.py";
          char *parallel_script = "lofreq2_call_pparallel.py";
          char *vcfset_script = "lofreq2_vcfset.py";
          char *vcfplot_script = "lofreq2_vcfplot.py";
          char *cluster_script = "lofreq2_cluster.py";
          char *script_to_call;

          if (strcmp(argv[1], "somatic") == 0) {
               script_to_call = somatic_script;
#ifdef PARALLEL
          } else if (strcmp(argv[1], "call-parallel") == 0) {
               script_to_call = parallel_script;
#endif
          } else if (strcmp(argv[1], "vcfset") == 0) {
               script_to_call = vcfset_script;
          } else if (strcmp(argv[1], "vcfplot") == 0) {
               script_to_call = vcfplot_script;
          } else if (strcmp(argv[1], "cluster") == 0) {
               script_to_call = cluster_script;
          } else {
               LOG_FATAL("%s\n", "Internal error: unknown option");
               return 1;
          }

          argv_execvp[0] = argv[0];
          for (i=2; i<argc; i++) {
               argv_execvp[i-1] = argv[i];
          }
          argv_execvp[i-1] = NULL; /* sentinel */
          if (execvp(script_to_call, argv_execvp)) {
               perror("Calling external LoFreq script via execvp failed");
               free(argv_execvp);
               return 1;
          } else {
               free(argv_execvp);
               return 0;
          }
#ifdef USE_ALNERRPROF
     } else if (strcmp(argv[1], "bamstats") == 0) {
          return main_bamstats(argc, argv);
#endif
     } else if (strcmp(argv[1], "plpsummary") == 0) {
          /* modify args to  main_call() */
          char **argv_tmp = calloc(argc+1, sizeof(char*));
          int i, rc;
          char plp_summary_arg[] = "--plp-summary-only";
          /*char bam_stats_arg[] = "--bam-stats";*/
          char *call_arg = NULL;
     
          if (strcmp(argv[1], "plpsummary") == 0) {
               call_arg = plp_summary_arg;
          } else {
               LOG_FATAL("%s\n", "Internal error: unknown option");
               return 1;
          }

          LOG_VERBOSE("'%s' is just an alias for %s call %s"
                      "  (ignoring all the snv-call specific options)\n", 
                      argv[1], BASENAME(argv[0]), call_arg);
          argv_tmp[0] = argv[0];
          argv_tmp[1] = "call";
          argv_tmp[2] = call_arg;
          for (i=2; i<argc; i++) {
               argv_tmp[i+1] = argv[i];
          }
#ifdef TRACE
          for (i=0; i<argc+1; i++) {
               LOG_FIXME("argv[%d] = %s\n", i, argv_tmp[i]);
          }
          exit(1);
#endif

          rc = main_call(argc+1, argv_tmp);

          free(argv_tmp);
          return rc;

     } else if (strcmp(argv[1], "version") == 0) {
          fprintf(stdout, "version: %s\ncommit: %s\nbuild-date: %s\n",
                  PACKAGE_VERSION, GIT_VERSION, __DATE__);
          return 0;

     } else {
          LOG_FATAL("Unrecognized command '%s'\n", argv[1]);
          return 1;
     }
     return 0;
}
