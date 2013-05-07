/* -*- c-file-style: "k&r"; indent-tabs-mode: nil; -*-
 *
 *
 * FIXME missing license
 *
 */

/* NOTE: this is by no means a generic vcf parser, since many
 * functions depends on the properties/format of your variants. Here,
 * we only use whatever is needed inside LoFreq
 */

#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "log.h"
#include "utils.h"
#include "vcf.h"

#define MISSING_VAL_STR "."
#define MISSING_VAL_CHAR '.'

#define LINE_BUF_SIZE 1<<12


int vcf_var_filtered(const var_t *var)
{
     if (! var->filter) {
          return 0;
     } else if (0 == strcmp(var->filter, ".")) {
          return 0;
     } else if (strlen(var->filter)>=4 && 0 == strcmp(var->filter, "PASS")) {
          return 0;
     } else {
          return 1;
     }
}

/* value for key will be stored in value if not NULL. value will NULL
 * if not found. Otherwise its allocated here and caller must free. FIXME
 * shoddily written */
int
vcf_var_has_info_key(char **value, const var_t *var, const char *key) {
     const char field_delimiter[] = ";";
     char *token;
     char *info;
     char *info_ptr;

     if (value) {
          (*value) = NULL;
     }

     if (! var->info) {
          return 0;
     } 

     info = strdup(var->info); /* strsep modifies string */
     info_ptr = info;
     while (NULL != (token = strsep(&info_ptr, field_delimiter))) {
          if (0 == strncasecmp(key, token, strlen(key))) {
               if (value) {
                    char *s = strchr(token, '=');
                    if (NULL != s) {
                         (*value) = strdup(s+1);
                    }
               }
               free(info);
               return 1;
          }
     }
     free(info);

     return 0;
}


void vcf_new_var(var_t **var)
{
     (*var) = malloc(sizeof(var_t));
     (*var)->chrom = NULL;
     (*var)->pos = -1;
     (*var)->id = NULL;
     (*var)->ref = MISSING_VAL_CHAR;
     (*var)->alt = MISSING_VAL_CHAR;
     (*var)->qual = -1;
     (*var)->filter = NULL;
     (*var)->info = NULL;

     (*var)->format = NULL;
     (*var)->num_samples = 0;
     (*var)->samples = NULL;
}


void vcf_free_var(var_t **var)
{
     int i;
     free((*var)->chrom);
     free((*var)->id);
     free((*var)->filter);
     free((*var)->info);

     free((*var)->format);
     for (i=0; i<(*var)->num_samples; i++) {
          free((*var)->samples[i]);
     }
     free((*var)->samples);

     free(*var);
}


void vcf_write_var(FILE *stream, const var_t *var)
{
     /* in theory all values are optional */

     fprintf(stream, "%s\t%ld\t%s\t%c\t%c\t",
             NULL == var->chrom ? MISSING_VAL_STR : var->chrom,
             var->pos + 1,
             NULL == var->id ? MISSING_VAL_STR : var->id,
             var->ref,
             var->alt);
     if (var->qual>-1) {
          fprintf(stream, "%d", var->qual);
     } else {
          fprintf(stream, "%c", MISSING_VAL_CHAR);
     }

     fprintf(stream, "\t%s\t%s\n",
             var->filter ? var->filter : MISSING_VAL_STR,
             var->info ? var->info : MISSING_VAL_STR);    

     /* FIXME format and samples not supported */
}


/* var->info allocated here. caller has to free */
void vcf_var_sprintf_info(var_t *var,
                         const int *dp, const float *af, const int *sb,
                         const dp4_counts_t *dp4,
                         const int indel, const int consvar)
{
     char buf[LINE_BUF_SIZE];          
     snprintf(buf, sizeof(buf)-32, /* leave some for INDEL and other flags below */
              "DP=%d;AF=%f;SB=%d;DP4=%d,%d,%d,%d",
              *dp, *af, *sb, dp4->ref_fw, dp4->ref_rv, dp4->alt_fw, dp4->alt_rv);
     if (indel) {
          sprintf(buf, "%s;INDEL", buf);
     }
     if (consvar) {
          sprintf(buf, "%s;CONSVAR", buf);
     }
     
     var->info = strdup(buf);

     /* FIXME format and samples not supported */
}


/* src can either be the program or the command. that's at least what
 * the vcftools folks do as well.
 */
void vcf_write_header(FILE *stream, const char *src, const char *reffa)
{
     char tbuf[9];
     struct tm tm;
     time_t t;

     t = time(0);
     localtime_r(&t, &tm);
     strftime(tbuf, 9, "%Y%m%d", &tm);
     
     fprintf(stream, "##fileformat=VCFv4.0\n");
     fprintf(stream, "##fileDate=%s\n", tbuf);
     if (src) {
          fprintf(stream, "##source=%s\n", src);
     }
     if (reffa) {
          fprintf(stream, "##reference=%s\n", reffa);
     }
     fprintf(stream, "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Raw Depth\">\n");
     fprintf(stream, "##INFO=<ID=AF,Number=1,Type=Float,Description=\"Allele Frequency\">\n");    
     fprintf(stream, "##INFO=<ID=SB,Number=1,Type=Integer,Description=\"Phred-scaled strand bias at this position\">\n");    

     fprintf(stream, "##INFO=<ID=DP4,Number=4,Type=Integer,Description=\"Counts for ref-forward bases, ref-reverse, alt-forward and alt-reverse bases\">\n");

     fprintf(stream, "##INFO=<ID=INDEL,Number=0,Type=Flag,Description=\"Indicates that the variant is an INDEL.\">\n");
     fprintf(stream, "##INFO=<ID=CONSVAR,Number=0,Type=Flag,Description=\"Indicates that the variant is a consensus variant (as opposed to a low frequency variant).\">\n");

     fprintf(stream, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");
}



/* parse header, i.e. meta info until and including header from vcf
 * file. will allocate memory for header. caller has to free. returns
 * 0 on success.
 */
int vcf_parse_header(char **header, FILE *stream)
{
     const char *HEADER_LINE = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO";
     char line[LINE_BUF_SIZE];

     /* make sure strlen below will work on header */
     (*header) = malloc(sizeof(char)); 
     (*header)[0] = '\0';

     while (NULL != fgets(line, sizeof(line), stream)) {
          (*header) = realloc((*header), (strlen(*header) + strlen(line) + 1 /* '\0' */) * sizeof(char));
          (void) strcat((*header), line);
          if (strlen(line) >= strlen(HEADER_LINE)) {
               if (0 == strncmp(line, HEADER_LINE, strlen(HEADER_LINE))) {
                    return 0;
               }
          }
     }
     LOG_WARN("%s\n", "Missing header line in vcf file.");
     return -1;
}


/* parse one variant from stream. returns +1 on EOF and -1 on error
 */
int parse_var(FILE *stream, var_t *var)
{
     const char delimiter[] = "\t";
     char *token;
     char line[LINE_BUF_SIZE];
     char *line_ptr;
     int field_no = 0;
     
     if (NULL == fgets(line, sizeof(line), stream)) {
          return 1;          
     }
     chomp(line);
     line_ptr = line;

     LOG_DEBUG("parsing line: %s\n", line);

     /* note: strsep modifies line_ptr */
     while (NULL != (token = strsep(&line_ptr, delimiter))) {
          field_no+=1;
          if (1 == field_no) {
               var->chrom = strdup(token);

          } else if (2 == field_no) {
               var->pos = atol(token)-1;               
          } else if (3 == field_no) {
               var->id = strdup(token);                      
               
          } else if (4 == field_no) {
               if (strlen(token)>1) {
                    LOG_FATAL("%s\n", "Only supporting one reference base in vcf");
               }
               var->ref = token[0];
               
          } else if (5 == field_no) {
               if (strlen(token)>1) {
                    LOG_FATAL("%s\n", "Only supporting one alt base in vcf");
               }
               var->alt = token[0];               
               
          } else if (6 == field_no) {
               var->qual = atoi(token);               

          } else if (7 == field_no) {
               var->filter = strdup(token);               
               
          } else if (8 == field_no) {
               var->info = strdup(token);               

#if 0
          } else if (9 == field_no) {
               var->format = strdup(token);               

          } else if (9 < field_no) {
               /* allocate mem for samples first */
               var->samples[field_no-10] = strdup(token);               
#else
          } else {
               LOG_WARN("%s\n", "Genotyping info in vcf not supported");
          }
#endif
     }

     if (field_no<8) {
          LOG_WARN("Parsing of variant incomplete. Only got %d fields. Need at least 8\n",
                   field_no);
          return -1;
     }

     return 0;
}


/* parse all variants from stream and return number of parsed vars or
 * -1 on error. memory for vars will be allocated here.
 */
int vcf_parse_vars(FILE *stream, var_t ***vars)
{
     int rc;
     int num_vars = 0;

     (*vars) = malloc(1 * sizeof(var_t*));

     while (! feof(stream)) {
          var_t *var;
          vcf_new_var(&var);
          rc = parse_var(stream, var);
          if (-1 == rc) {
               int i;
               LOG_FATAL("%s\n", "Parsing error");               
               for (i=0; i<num_vars; i++) {
                    vcf_free_var(&(*vars)[i]);
               }
               vcf_free_var(&var);
               free((*vars));
               return -1;
          }
          if (1 == rc) {/* EOF */
               free(var);
               break;
          }
          num_vars += 1;
          (*vars) = realloc((*vars), num_vars * sizeof(var_t*));
          (*vars)[num_vars-1] = var;
#if 0
          LOG_DEBUG("(*vars)[num_vars-1 = %d] = \n", num_vars-1);
          vcf_write_var(stderr, (*vars)[num_vars-1]);
#endif
     }

     return num_vars;
}


#ifdef VCF_MAIN


/* 
gcc -pedantic -Wall -g -std=gnu99 -O2 -DVCF_MAIN -o vcf_main vcf.c utils.c log.c 
valgrind --tool=memcheck --leak-check=full --show-reachable=yes --track-origins=yes ./vcf_main example.vcf
*/
int main(int argc, char *argv[]) {
     char *header;
     var_t **vars = NULL;
     FILE *fh;
     int num_vars = 0;
     int i;

#if 0
     debug = 1;
     verbose = 1;
#endif


     if (argc < 2) {
          LOG_FATAL("%s\n", "Need vcf file as input");
          return 1;
     }


     fh = fopen(argv[1], "r");
     if (0 !=  vcf_parse_header(&header, fh)) {
          LOG_FATAL("%s\n", "vcf_parse_header() failed");
          free(header);
          return 1;
     } 
     fprintf(stdout, "HEADER start:\n");
     vcf_write_header(stdout, NULL, NULL);
     fprintf(stdout, "HEADER END\n");

     free(header);
     


     if (-1 == (num_vars = vcf_parse_vars(fh, &vars))) {
          LOG_FATAL("%s\n", "vcf_parse_vars() failed");
          return 1;
     }
     fprintf(stdout, "\nWRITING %d vars\n", num_vars);
     for (i=0; i<num_vars; i++) {
          vcf_write_var(stdout, vars[i]);
     }

     for (i=0; i<num_vars; i++) {
          vcf_free_var(& vars[i]);
     }
     free(vars);


     fclose(fh);
     LOG_VERBOSE("%s\n", "successful exit");

     return 0;
}
#endif
