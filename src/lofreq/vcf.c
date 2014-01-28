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
#include "uthash.h"

#include "log.h"
#include "utils.h"
#include "vcf.h"


#define LINE_BUF_SIZE 1<<12

/* this is the actual header. all the other stuff is actually called meta-info */
const char *HEADER_LINE = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";


int warned_one_alt_base_support = 0;
int warned_one_ref_base_support = 0;

void 
var_hash_free_table(var_hash_t *var_hash)
{
     var_hash_t *cur, *tmp;
     if (NULL == var_hash) {
          return;
     }
     HASH_ITER(hh, var_hash, cur, tmp) {
#ifdef TRACE
          LOG_ERROR("Freeing %s\n", cur->key);
#endif
          HASH_DEL(var_hash, cur);
          var_hash_free_elem(cur);
     }
}


void
var_hash_free_elem(var_hash_t *hash_elem_ptr) 
{
     vcf_free_var(& hash_elem_ptr->var);
     free(hash_elem_ptr->key);
     free(hash_elem_ptr);
}

/* key and var will not be copied ! */
void
var_hash_add(var_hash_t **var_hash, char *key, var_t *var)
{
    var_hash_t *vh_elem = NULL;
    vh_elem = (var_hash_t *) malloc(sizeof(var_hash_t));
    vh_elem->key = key;
    vh_elem->var = var;

    HASH_ADD_KEYPTR(hh, (*var_hash), vh_elem->key, strlen(vh_elem->key), vh_elem);
}



/* key is allocated here and has to be freed by called */
void
vcf_var_key(char **key, var_t *var)
{
     int bufsize = strlen(var->chrom)+16;
     (*key) = malloc(bufsize *sizeof(char));
     snprintf(*key, bufsize, "%s %ld %c %c", var->chrom, var->pos+1, 
              var->ref, var->alt);
     /* pos+1 make terminal output easier */
}

/* as above but only using chrom and pos */
void
vcf_var_key_simple(char **key, var_t *var)
{
     int bufsize = strlen(var->chrom)+16;
     (*key) = malloc(bufsize *sizeof(char));
     snprintf(*key, bufsize, "%s %ld", var->chrom, var->pos+1);
     /* pos+1 make terminal output easier */     
}



int
vcf_file_seek(vcf_file_t *f, long int offset, int whence) 
{
     if (f->gz) {
          return gzseek(f->fh_gz, offset, whence);
     } else {
          return fseek(f->fh, offset, whence);
     }
}


/* returns 0 on success. non-zero otherwise */
int
vcf_file_open(vcf_file_t *f, const char *path, const int gzip, char mode) 
{
     /* 1. gzip in/out should be mode b, right? if not, the code can be
      * simplified 
      *
      * 2. FIXME handle stdin and stdout
      */

     if (mode!='r' && mode!='w') {
          LOG_FATAL("Internal error: unknown mode %c\n", mode);
          return -1;
     }

     if (path[0] != '-' && mode=='r') {
          if (! file_exists(path) || is_dir(path)) {
               LOG_ERROR("VCF file %s does not exist\n", path);
               return -1;
          }
     }

     if (gzip) {
          if (path[0] == '-') {
               LOG_FIXME("%s\n", "gzip support for stdin/stdout not implemented yet");
               return -1;
          }
          f->gz = 1;
          f->fh = NULL;
          if (mode=='r') {
               f->fh_gz = gzopen(path, "rb");
          } else if (mode=='w') {
               f->fh_gz = gzopen(path, "wb");
          }

     } else {
          f->gz = 0;
          f->fh_gz = NULL;
          if (mode=='r') {
               if (path[0] == '-') {
                    f->fh = stdin;
               } else {
                    f->fh = fopen(path, "r");
               }
          } else if (mode=='w') {
               if (path[0] == '-') {
                    f->fh = stdout;
               } else {
                    f->fh = fopen(path, "w");
               }
          }
     }     

     if (! f->fh && ! f->fh_gz) {
          return -1;
     } else {
          return 0;
     }
}


int
vcf_file_close(vcf_file_t *f) 
{
     /* FIXME handle stdin and stdout */

     if (f->gz) {          
          return gzclose(f->fh_gz);
     } else {
          if (f->fh!=stdout) {
               return fclose(f->fh);
          } else {
               return 0;
          }
     }
}

char *
vcf_file_gets(vcf_file_t *f, int len, char *line) 
{
     if (f->gz) {
          return gzgets(f->fh_gz, line, len);

     } else {
          return fgets(line, len, f->fh);
     }
}



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
 * if not found. Otherwise its allocated here and caller must free.
 * FIXME shoddily written and we should use a hash for info key:val
 * pairs anyway */
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
     (*var)->ref = VCF_MISSING_VAL_CHAR;
     (*var)->alt = VCF_MISSING_VAL_CHAR;
     (*var)->qual = -1; /* -1 == missing */
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


void vcf_write_var(vcf_file_t *vcf_file, const var_t *var)
{
     /* in theory all values are optional */

     VCF_PRINTF(vcf_file, "%s\t%ld\t%s\t%c\t%c\t",
             NULL == var->chrom ? VCF_MISSING_VAL_STR : var->chrom,
             var->pos + 1,
             NULL == var->id ? VCF_MISSING_VAL_STR : var->id,
             var->ref,
             var->alt);
     if (var->qual>-1) {
          VCF_PRINTF(vcf_file, "%d", var->qual);
     } else {
          VCF_PRINTF(vcf_file, "%c", VCF_MISSING_VAL_CHAR);
     }

     VCF_PRINTF(vcf_file, "\t%s\t%s",
             var->filter ? var->filter : VCF_MISSING_VAL_STR,
             var->info ? var->info : VCF_MISSING_VAL_STR);

     if (var->format) {
          int i=0;
          VCF_PRINTF(vcf_file, "\t%s", var->format);
          for (i=0; i<var->num_samples; i++) {
               VCF_PRINTF(vcf_file, "\t%s", var->samples[i]);
          }
     }
     VCF_PRINTF(vcf_file, "\n");
}


void vcf_var_add_to_info(var_t *var, const char *info_str)
{
     var->info = realloc(var->info,
                         (strlen(var->info) + strlen(info_str)
                          + 1/*;*/ + 1/*\0*/) * sizeof(char));
     if (strlen(var->info)) {
          (void) strcat(var->info, ";");
     }
     (void) strcat(var->info, info_str);
}


void vcf_var_add_to_filter(var_t *var, const char *filter_name)
{
     if (var->filter) {
          /* clear field, if PASSED or missing  */
          if ((strlen(var->filter)>=4 && 0 == strcmp(var->filter, "PASS"))
              ||
              (strlen(var->filter) && var->filter[0] == VCF_MISSING_VAL_CHAR)) {
               free(var->filter);
               var->filter = NULL;
          }
     }

     if (! var->filter) {/* could have been freed above so don't else if */
          var->filter = malloc(1 * sizeof(char));
          var->filter[0] = '\0';
     }

     /* realloc */
     if (var->filter) {
          var->filter = realloc(var->filter,
                                (strlen(var->filter) + strlen(filter_name)
                                + 1/*;*/ + 1/*\0*/) * sizeof(char));
     }

     if (var->filter==NULL) {
          fprintf(stderr, "FATAL: couldn't allocate memory at %s:%s():%d\n",
                  __FILE__, __FUNCTION__, __LINE__);
     }

     /* add */
     if (strlen(var->filter)) {
          (void) strcat(var->filter, ";");
     }
     (void) strcat(var->filter, filter_name);
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


void vcf_write_header(vcf_file_t *vcf_file, const char *header)
{
     VCF_PRINTF(vcf_file, "%s", header);
}

/* src can either be the program or the command. that's at least what
 * the vcftools folks do as well.
 */
void vcf_write_new_header(vcf_file_t *vcf_file, const char *src, const char *reffa)
{
     char tbuf[9];
     struct tm tm;
     time_t t;

     t = time(0);
     localtime_r(&t, &tm);
     strftime(tbuf, 9, "%Y%m%d", &tm);

     VCF_PRINTF(vcf_file, "##fileformat=VCFv4.0\n");
     VCF_PRINTF(vcf_file, "##fileDate=%s\n", tbuf);
     if (src) {
          VCF_PRINTF(vcf_file, "##source=%s\n", src);
     }
     if (reffa) {
          VCF_PRINTF(vcf_file, "##reference=%s\n", reffa);
     }
     VCF_PRINTF(vcf_file, "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Raw Depth\">\n");
     VCF_PRINTF(vcf_file, "##INFO=<ID=AF,Number=1,Type=Float,Description=\"Allele Frequency\">\n");
     VCF_PRINTF(vcf_file, "##INFO=<ID=SB,Number=1,Type=Integer,Description=\"Phred-scaled strand bias at this position\">\n");
     VCF_PRINTF(vcf_file, "##INFO=<ID=DP4,Number=4,Type=Integer,Description=\"Counts for ref-forward bases, ref-reverse, alt-forward and alt-reverse bases\">\n");
     VCF_PRINTF(vcf_file, "##INFO=<ID=INDEL,Number=0,Type=Flag,Description=\"Indicates that the variant is an INDEL.\">\n");
     VCF_PRINTF(vcf_file, "##INFO=<ID=CONSVAR,Number=0,Type=Flag,Description=\"Indicates that the variant is a consensus variant (as opposed to a low frequency variant).\">\n");
     VCF_PRINTF(vcf_file, "%s", HEADER_LINE);
}



/* parse header, i.e. meta info until and including header from vcf
 * file. will allocate memory for header. caller has to free. returns
 * 0 on success.
 */
int vcf_parse_header(char **header, vcf_file_t *vcf_file)
{
     char line[LINE_BUF_SIZE];
     
     /* make sure strlen below will work on header */
     (*header) = malloc(sizeof(char));
     (*header)[0] = '\0';

     while (NULL != vcf_file_gets(vcf_file, sizeof(line), line)) {

#ifdef TRACE
          fprintf(stderr, "Got line %s\n", line);
#endif
          (*header) = realloc((*header), (strlen(*header) + strlen(line) + 1 /* '\0' */) * sizeof(char));
          (void) strcat((*header), line);
          if (strlen(line) >= strlen(HEADER_LINE)-1) {
               if (0 == strncmp(line, HEADER_LINE, strlen(HEADER_LINE)-1)) {
                    return 0;
               }
          }
     }

     free(*header);
     LOG_WARN("%s\n", "Missing header line in vcf file.");
     return -1;
}


int vcf_skip_header(vcf_file_t *vcf_file)
{
     char *vcf_header;
     if (0 !=  vcf_parse_header(&vcf_header, vcf_file)) {
          if (vcf_file_seek(vcf_file, 0, SEEK_SET)) {
               LOG_FATAL("%s\n", "Couldn't rewind file to parse variants"
                        " after header parsing failed");
              return -1;
         }
     } else {
          free(vcf_header);
     }
     return 0;
}



/* parse one variant from stream. returns +1 on EOF and -1 on error
 */
int vcf_parse_var(vcf_file_t *vcf_file, var_t *var)
{
     const char delimiter[] = "\t";
     char *token;
     char line[LINE_BUF_SIZE];
     char *line_ptr;
     int field_no = 0;

     if (NULL == vcf_file_gets(vcf_file, sizeof(line), line)) {
          return 1;
     }
     chomp(line);
     line_ptr = line;

#if 0
     LOG_DEBUG("parsing line: %s\n", line);
#endif

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
               if (! warned_one_ref_base_support && strlen(token)>1) {
                    LOG_WARN("%s\n", "Only supporting one reference base in vcf");
                    warned_one_ref_base_support = 1;
               }
               var->ref = token[0];

          } else if (5 == field_no) {
               if (! warned_one_alt_base_support && strlen(token)>1) {
                    LOG_WARN("%s\n", "Only supporting one alt base in vcf");
                    warned_one_alt_base_support = 1;
               }
               var->alt = token[0];

          } else if (6 == field_no) {
               if (token[0]==VCF_MISSING_VAL_CHAR) {
                    var->qual = -1;
               } else {
                    var->qual = atoi(token);
               }

          } else if (7 == field_no) {
               var->filter = strdup(token);

          } else if (8 == field_no) {
               var->info = strdup(token);

          } else if (9 == field_no) {
               var->format = strdup(token);

          } else if (field_no > 9) {
               assert(field_no-10 == var->num_samples);
               var->num_samples += 1;
               var->samples = realloc(var->samples, var->num_samples * sizeof(char*));
               var->samples[var->num_samples-1] = strdup(token);
          }
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
int vcf_parse_vars(var_t ***vars, vcf_file_t *vcf_file, int only_passed)
{
     int rc;
     int num_vars = 0;

     (*vars) = malloc(1 * sizeof(var_t*));

     while (! VCF_EOF(vcf_file)) { 
          var_t *var;
          vcf_new_var(&var);
          rc = vcf_parse_var(vcf_file, var);
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

          if (only_passed==1) {
               if (vcf_var_filtered(var)) {
                    vcf_free_var(&var);
                    continue;
               }
          }
          num_vars += 1;
          (*vars) = realloc((*vars), num_vars * sizeof(var_t*));
          (*vars)[num_vars-1] = var;
          if (verbose && num_vars && num_vars%1000000==0) {
               LOG_VERBOSE("Still alive and happily parsing var %d\n", num_vars);
          }
#if 0
          LOG_DEBUG("(*vars)[num_vars-1 = %d] = \n", num_vars-1);
          vcf_write_var(stderr, (*vars)[num_vars-1]);
#endif
     }

     return num_vars;
}



/* info needs to be terminated with a newline character */
void vcf_header_add(char **header, const char *info)
{
     char *token;
     int pos;

     /* make sure to insert before HEADER_LINE */

     token = strstr(*header, HEADER_LINE);
     if (! token) {
          LOG_WARN("%s\n", "Can't add info to empty header");
          return;
     }
     pos = (int)(token - (*header));

     *header = realloc(*header, (strlen(*header) + strlen(info) + 1) * sizeof(char));

#if 0
     LOG_FIXME("header-len=%d; info len=%d; alloc=%d; pos=%d\n",
        strlen(*header), strlen(info), (strlen(*header) + strlen(info) + 1), pos);
#endif

     (*header)[pos] = '\0'; /* can't just: token[0] = '\0'; since that would work on a copy?! */
     (void) strcat(*header, info);
     (void) strcat(*header, HEADER_LINE);
     return;
}



#ifdef VCF_MAIN


/*
gcc  -Wall -g -std=gnu99 -O2 -DVCF_MAIN -o vcf_main vcf.c utils.c log.c -lz
valgrind --tool=memcheck --leak-check=full --show-reachable=yes --track-origins=yes ./vcf_main example.vcf
*/
int main(int argc, char *argv[]) {
     char *header;
     var_t **vars = NULL;
     vcf_file_t vcf_file_in, vcf_file_out;
     int num_vars = 0;
     int i;
     char *path_in, *path_out;
     int gzip_in, gzip_out = 0;
#if 0
     debug = 1;
     verbose = 1;
#endif


     if (argc < 3) {
          LOG_FATAL("%s\n", "Need two args: vcf-in vcf-out");
          return 1;
     }
     path_in = argv[1];
     path_out = argv[2];
     
     if (HAS_GZIP_EXT(path_in)) {
          gzip_in = 1;
     } else {
          gzip_in = 0;
     }
     LOG_FIXME("Using %s (%s gzipped)\n", path_in, gzip_in ? "is" : "not");
     if (vcf_file_open(& vcf_file_in, path_in, gzip_in, 'r')) {
          LOG_FATAL("%s\n", "vcf_file_open() failed");
          exit(1);
     }

     if (HAS_GZIP_EXT(path_out)) {
          gzip_out = 1;
     } else {
          gzip_out = 0;
     }
     LOG_FIXME("Using %s (%s gzipped)\n", path_out, gzip_out ? "is" : "not");
     if (vcf_file_open(& vcf_file_out, path_out, gzip_out, 'w')) {
          LOG_FATAL("%s\n", "vcf_file_open() failed");
          exit(1);
     }


     if (0 !=  vcf_parse_header(&header, & vcf_file_in)) {
          LOG_FATAL("%s\n", "vcf_parse_header() failed");
          free(header);
          return 1;
     }
     vcf_write_new_header(& vcf_file_out, NULL, NULL);
     free(header);


     if (-1 == (num_vars = vcf_parse_vars(& vcf_file_in, &vars))) {
          LOG_FATAL("%s\n", "vcf_parse_vars() failed");
          return 1;
     }
     fprintf(stdout, "Wrote %d vars to output\n", num_vars);
     for (i=0; i<num_vars; i++) {
          vcf_write_var(& vcf_file_out, vars[i]);
     }

     for (i=0; i<num_vars; i++) {
          vcf_free_var(& vars[i]);
     }
     free(vars);

     vcf_file_close(& vcf_file_in);
     vcf_file_close(& vcf_file_out);

     LOG_VERBOSE("%s\n", "successful exit");

     return 0;
}
#endif
