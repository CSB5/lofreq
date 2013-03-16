/* -*- c-file-style: "k&r" -*-
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "log.h"
#include "utils.h"
#include "lofreq_snpcaller.h"


static void usage(const char *myname)
{
     fprintf(stderr, "%s: Fast and sensitive inference of single-nucleotide variants\n", PACKAGE_NAME);
     /*fprintf(stderr, "Version %s\n", PACKAGE_VERSION);*/
     fprintf(stderr, "\n");
     fprintf(stderr, "Usage: %s <command> [options], where command is one of:\n", myname);
     fprintf(stderr, "  call        : call variants\n");
     fprintf(stderr, "  filter      : filter variants\n");
     fprintf(stderr, "  plp_summary : print pileup summary\n");
     fprintf(stderr, "  version     : prints version\n");
     fprintf(stderr, "\n");
}

int main(int argc, char *argv[])
{
     if (argc < 2) {
          usage(BASENAME(argv[0]));
          return 1;
     }
     if (strcmp(argv[1], "call") == 0)  {
          return main_call(argc-1, argv+1);

     } else if (strcmp(argv[1], "filter") == 0) {
          char **argv_execvp = calloc(argc-1, sizeof(char*));
          int i;
          argv_execvp[0] = argv[0];
          for (i=2; i<argc; i++) {
               argv_execvp[i-1] = argv[i];
          }
          if (execvp("lofreq2_filter.py", argv_execvp)) {
               perror("Calling lofreq2_filter.py via execvp failed");
               free(argv_execvp);
               return -1;
          } else {
               free(argv_execvp);
               return 0;
          }

     } else if (strcmp(argv[1], "plp_summary") == 0) {
          /* call self with command 'call' but add --plp_summary */
          char **argv_tmp = calloc(argc, sizeof(char*));
          int i, rc;
          fprintf(stderr, "NOTE: the plp_summary command is just an alias for %s call --plp-summary-only (ignoring all the snv-call specific options)\n", BASENAME(argv[0]));
          argv_tmp[0] = strdup("call");
          for (i=2; i<argc; i++) {
               argv_tmp[i-1] = argv[i];
          }
          argv_tmp[argc-1] = "--plp-summary-only";
          rc = main_call(argc, argv_tmp);
          free(argv_tmp);
          return rc;

     } else if (strcmp(argv[1], "version") == 0) {
          fprintf(stdout, "%s\n", PACKAGE_VERSION);
          return 0;
     } else {
          LOG_FATAL("Unrecognized command '%s'\n", argv[1]);
          return 1;
     }
     return 0;
}
