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

#include "log.h"
#include "utils.h"
#include "lofreq_snpcaller.h"
#include "lofreq_uniq.h"




/* prepend dirname(argv0) and python source dir to PATH making sure
 * that package works even without properly installing it
 */
void 
add_local_dir_to_path(char *argv0) {
     const char *PATH_NAME = "PATH";
     char *path_var = NULL;
     char *old_path = NULL;
     char *argv0_cp = NULL;
     char *dirname_argv0 = NULL;
     char lofreq_script_rel[] = "../lofreq_python/scripts/lofreq2_filter.py";
     char *lofreq_script_abs = NULL;

     argv0_cp = strdup(argv0); /* dirname might change contents */
     if (NULL == (dirname_argv0 = strdup(dirname(argv0_cp)))) {
          free(argv0_cp);
          return;
     }

     if (NULL == (lofreq_script_abs = strdup(dirname_argv0))) {
          free(dirname_argv0);
          free(argv0_cp);
          return;
     }

     if (NULL == join_paths(&lofreq_script_abs, lofreq_script_rel)) {
#if 0
          LOG_WARN("join_paths %s and %s failed\n", lofreq_script_abs, lofreq_script_rel);
#endif
          free(lofreq_script_abs);
          free(dirname_argv0);
          free(argv0_cp);
          return;          
     }
     /* check also done by realpath in join_path but doesn't hurt to
      * check again */
     if (! file_exists(lofreq_script_abs)) {
#if 0
          LOG_WARN("%s doesnt' exist\n", lofreq_script_abs);
#endif
          free(lofreq_script_abs);
          free(dirname_argv0);
          free(argv0_cp);
          return;
     }

     path_var = strdup(dirname(lofreq_script_abs));
#if 0
     LOG_WARN("Adding local source directory %s to PATH\n", path_var);
#endif

     old_path = getenv(PATH_NAME);
     if (NULL == old_path) {
          setenv(PATH_NAME, path_var, 1);
     } else {
          path_var = realloc(
               path_var, (strlen(path_var) + 1 + strlen(old_path) + 1) * sizeof(char));
          (void) strcat(path_var, ":");
          (void) strcat(path_var, old_path);
          setenv(PATH_NAME, path_var, 1);
     }


#ifdef FIXME_NO_NEED
     free(path_var);
     old_path = getenv(PATH_NAME);
     path_var = strdup(dirname(argv0));
     LOG_VERBOSE("Adding argv0 directory %s to PATH\n", path_var);
     path_var = realloc(
          path_var, (strlen(path_var) + 1 + strlen(old_path) + 1) * sizeof(char));
     (void) strcat(path_var, ":");
     (void) strcat(path_var, old_path);
     setenv(PATH_NAME, path_var, 1);
     LOG_DEBUG("New PATH = %s\n", path_var);
#endif

     free(path_var);
     free(lofreq_script_abs);
     free(dirname_argv0);
     free(argv0_cp);
}


static void usage(const char *myname)
{
     fprintf(stderr, "\n%s: Fast and sensitive inference of single-nucleotide variants\n", PACKAGE_NAME);
     /*fprintf(stderr, "Version %s\n", PACKAGE_VERSION);*/
     fprintf(stderr, "\n");
     fprintf(stderr, "Usage: %s <command> [options]\n\n", myname);
     fprintf(stderr, "  main commands:\n");
     fprintf(stderr, "    call        : call variants\n");
     fprintf(stderr, "    somatic     : call somatic variants\n\n");
     fprintf(stderr, "  other commands:\n");
     fprintf(stderr, "    filter      : filter variants\n");
     fprintf(stderr, "    peek        : check properties of BAM file\n");
     fprintf(stderr, "    uniq        : test whether SNVs predicted in one sample couldn't be predicted in other\n");
     fprintf(stderr, "    plp_summary : print pileup summary per position\n");
     fprintf(stderr, "    vcfset      : VCF set operations\n");
     fprintf(stderr, "    version     : print version string\n");
     fprintf(stderr, "\n");
}


int main(int argc, char *argv[])
{
    /* FIXME warn if link
     * #include <sys/stat.h>
     * int x = lstat("junklink", &buf);
     * if if (S_ISLNK(buf.st_mode))
     * warning things might not work as expected
     * else 
     */
     add_local_dir_to_path(argv[0]);

     if (argc < 2) {
          usage(BASENAME(argv[0]));
          return 1;
     }
     if (strcmp(argv[1], "call") == 0)  {
          return main_call(argc, argv);

     } else if (strcmp(argv[1], "uniq") == 0)  {
          return main_uniq(argc, argv);

     } else if (strcmp(argv[1], "peek") == 0 ||
                strcmp(argv[1], "check") == 0 ||
                strcmp(argv[1], "inspect") == 0 ||
                strcmp(argv[1], "doctor") == 0 ||
                strcmp(argv[1], "run-me-first") == 0)  {
          LOG_FIXME("%s\n", "NOT IMPLEMENT YET: sq_list, has_delqs, has_insqs, paired_reads...\n");
          return -1;


     } else if (strcmp(argv[1], "filter") == 0 || 
                strcmp(argv[1], "somatic") == 0 ||
                strcmp(argv[1], "vcfset") == 0) {
          char **argv_execvp = calloc(argc, sizeof(char*));
          int i;
          char *filter_script = "lofreq2_filter.py";
          char *somatic_script = "lofreq2_somatic.py";
          char *vcfset_script = "lofreq2_vcfset.py";
          char *script_to_call;;

          if (strcmp(argv[1], "filter") == 0) {
               script_to_call = filter_script;
          } else if (strcmp(argv[1], "somatic") == 0) {
               script_to_call = somatic_script;
          } else if (strcmp(argv[1], "vcfset") == 0) {
               script_to_call = vcfset_script;
          } else {
               LOG_FATAL("%s\n", "Internal error: unknown option");
               return -1;
          }

          argv_execvp[0] = argv[0];
          for (i=2; i<argc; i++) {
               argv_execvp[i-1] = argv[i];
          }
          argv_execvp[i-1] = NULL; /* sentinel */
          if (execvp(script_to_call, argv_execvp)) {
               perror("Calling external LoFreq script via execvp failed");
               free(argv_execvp);
               return -1;
          } else {
               free(argv_execvp);
               return 0;
          }

     } else if (strcmp(argv[1], "plp_summary") == 0) {
          /* use main_call() but add --plp_summary */
          char **argv_tmp = calloc(argc+1, sizeof(char*));
          int i, rc;
          LOG_VERBOSE("'plp_summary' is just an alias for %s call --plp-summary-only"
                      "  (ignoring all the snv-call specific options)\n", BASENAME(argv[0]));
          argv_tmp[0] = argv[0];
          argv_tmp[1] = "call";
          argv_tmp[2] = "--plp-summary-only";
          for (i=2; i<argc; i++) {
               argv_tmp[i+1] = argv[i];
          }
#if 0
          for (i=0; i<argc+1; i++) {
               LOG_FIXME("New arg %d: %s\n", i, argv_tmp[i]);
          }
#endif
          rc = main_call(argc+1, argv_tmp);
          free(argv_tmp);
          return rc;

     } else if (strcmp(argv[1], "version") == 0) {
          fprintf(stdout, "version: %s\ncommit: %s\n", PACKAGE_VERSION, GIT_VERSION);
          return 0;

     } else {
          LOG_FATAL("Unrecognized command '%s'\n", argv[1]);
          return 1;
     }
     return 0;
}
