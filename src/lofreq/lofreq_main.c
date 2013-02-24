/* -*- c-file-style: "k&r" -*-
 *
 */
#include <stdio.h>
#include <string.h>

#include "lofreq_snpcaller.h"
#include "log.h"

#define BASENAME(x) strrchr((x), '/') ? strrchr((x), '/')+1 : (x)
static void usage(const char *myname)
{
     fprintf(stderr, "%s: Fast and sensitive inference of single-nucleotide variants\n", PACKAGE_NAME);
     fprintf(stderr, "Version %s\n", PACKAGE_VERSION);
     fprintf(stderr, "\n");
     fprintf(stderr, "Usage: %s <command> [options], where command is one of:\n", myname);
     fprintf(stderr, "  call   : call variants\n");
     fprintf(stderr, "  filter : filter variants\n");
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
          LOG_FIXME("%s\n", "not implemented yet");
          return -1;

     } else {
          LOG_FATAL("Unrecognized command '%s'\n", argv[1]);
          return 1;
     }
     return 0;
}
