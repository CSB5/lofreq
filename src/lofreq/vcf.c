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


/* NOTE: this is by no means a generic vcf parser, since many
 * functions depends on the properties/format of your variants. Here,
 * we only use whatever is needed inside LoFreq
 */

#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "htslib/bgzf.h"
#include "htslib/kstring.h"
#include "htslib/kseq.h"
#include "htslib/tbx.h"

#include "uthash.h"

#include "log.h"
#include "utils.h"
#include "vcf.h"
#include "defaults.h"

#define LINE_BUF_SIZE 1<<12


/* this is the actual header. all the other stuff is actually called meta-info 
 * note, newline character is missing here
 */
const char *VCF_HEADER = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO";



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


/* FIXME key and var will not be copied, i.e. don't free for now */
void
var_hash_add(var_hash_t **var_hash, char *key, var_t *var)
{
    var_hash_t *vh_elem = NULL;
    var_hash_t *match = NULL;

    HASH_FIND_STR((*var_hash), key, match);
    if (match) {
         LOG_DEBUG("Already got a variant match for key '%s'. Will keep the old one.\n", key);
         return;
    }

    vh_elem = (var_hash_t *) malloc(sizeof(var_hash_t));
    vh_elem->key = key;
    vh_elem->var = var;
    /* FIXME should we test for existance first? */
    HASH_ADD_KEYPTR(hh, (*var_hash), vh_elem->key, strlen(vh_elem->key), vh_elem);
}


/* key is allocated here and has to be freed by called */
void
vcf_var_key(char **key, var_t *var)
{
     int bufsize = strlen(var->chrom)+16;
     assert(var->ref && var->alt);
     (*key) = malloc(bufsize *sizeof(char));
     snprintf(*key, bufsize, "%s %ld %s %s", var->chrom, var->pos+1, 
              var->ref, var->alt);
     /* pos+1 make terminal output easier */
}

/* as above but only using chrom and pos */
void
vcf_var_key_pos_only(char **key, var_t *var)
{
     int bufsize = strlen(var->chrom)+16;
     (*key) = malloc(bufsize *sizeof(char));
     snprintf(*key, bufsize, "%s %ld", var->chrom, var->pos+1);
     /* pos+1 make terminal output easier */     
}


int vcf_printf(vcf_file_t *f, char *fmt, ...)
{
     /* sadly there is no gzvprintf */
     char buf[64000];/* needs to be able to hold header */
     va_list args;
     int len;

     va_start(args, fmt);
     len = vsnprintf(buf, 64000, fmt, args);    
     va_end(args);

     if (len>=64000) {
          LOG_WARN("%s\n", "Truncated vcf_printf");
     }
     if (f->is_bgz) {
          return bgzf_write(f->fh_bgz, buf, strlen(buf));
     } else {
          return fputs(buf, f->fh);
     }
}

int
vcf_file_seek(vcf_file_t *f, long int offset, int whence) 
{
     if (f->is_bgz) {
          return bgzf_seek(f->fh_bgz, offset, whence);
     } else {
          return fseek(f->fh, offset, whence);
     }
}


/* returns 0 on success. non-zero otherwise */
int
vcf_file_open(vcf_file_t *f, const char *path, const int bgzip, char mode) 
{
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

     f->path = strdup(path);
     f->mode =mode;
     
     if (bgzip) {
          if (path[0] == '-') {
               LOG_FIXME("%s\n", "bgzip support for stdin/stdout not implemented yet");
               return -1;
          }
          f->is_bgz = 1;
          f->fh = NULL;
          if (mode=='r') {
               f->fh_bgz = bgzf_open(path, "rb");
          } else if (mode=='w') {
               f->fh_bgz = bgzf_open(path, "wb");
          }

     } else {
          f->is_bgz = 0;
          f->fh_bgz = NULL;
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

     if (! f->fh && ! f->fh_bgz) {
          return -1;
     } else {
          return 0;
     }
}


int
vcf_file_flush(vcf_file_t *f)
{
     if (f->is_bgz) {          
          return bgzf_flush(f->fh_bgz);
     } else {
          return fflush(f->fh);
     }
}


/* note: tries to tabix index and also frees path */
int
vcf_file_close(vcf_file_t *f) 
{
     int rc = 0;
     if (f->is_bgz) {          
          rc = bgzf_close(f->fh_bgz);
          if (rc==0 && f->mode=='w' && f->path && f->path[0] != '-') {
               int min_shift = -1;
               tbx_conf_t conf = tbx_conf_vcf;
               rc = tbx_index_build(f->path, min_shift, &conf);
               if (rc) {
                    LOG_WARN("indexing of %s failed\n", f->path);
               }
          }
     } else {
          if (f->fh!=stdout) {
               rc = fclose(f->fh);
          } else {
               rc = 0;
          }
     }
     free(f->path);
     return rc;
}


/* returns NULL on error or EOF */
char *
vcf_file_gets(vcf_file_t *f, int len, char *line) 
{
     if (f->is_bgz) {
          kstring_t str = {0, 0, 0};
          if (bgzf_getline(f->fh_bgz, '\n', &str) > 0) {
               /* will get errors like
                  [E::get_intv] failed to parse TBX_VCF, was wrong -p [type] used?
                  The offending line was: "19,0,1"
                  on just gzipped data. not sure how to catch this. the following is a paranoia check
               */
               if (str.l<1) {
                    return NULL;
               }
               strncpy(line, str.s, len-2);
               /* behave like fgets and keep newline */
               line[strlen(line)] = '\n';
               line[strlen(line)+1] = '\0';
               
               free(str.s);
               return line;
          } else {
               return NULL;
          }

     } else {
          return fgets(line, len, f->fh);
     }
}


int vcf_var_filtered(const var_t *var)
{
     if (! var->filter) {
          return 0;
     } else if (0 == strcmp(var->filter, VCF_MISSING_VAL_STR)) {
          return 0;
     } else if (strlen(var->filter)>=4 && 0 == strcmp(var->filter, "PASS")) {
          return 0;
     } else {
          return 1;
     }
}

int vcf_var_is_indel(const var_t *var)
{
     if (strlen(var->ref)>1 ||
         strlen(var->alt)>1 ||
         vcf_var_has_info_key(NULL, var, "INDEL")) {
          return 1;
     } else {
          return 0;
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

     if (! var->info || ! key) {
          return 0;
     }
     if (strlen(var->info)<2) {
          return 0;
     }
     info = strdup(var->info);
     if (! info) {
          LOG_FATAL("%s\n", "insufficient memory");
          exit(1);
     }
     info_ptr = info;
     token = info;
     /* note: strsep modifies ptr. see also
      * http://stackoverflow.com/questions/21383082/parsing-a-string-in-c-with-strsep-alternative-methods */
     while (token) {
          strsep(&info_ptr, field_delimiter);
          /*fprintf(stderr, "token=%s key=%s\n", token, key);*/
          if (0 == strncasecmp(key, token, MIN(strlen(token), strlen(key)))) {
               if (value) {
                    char *s = strchr(token, '=');
                    if (NULL != s) {
                         (*value) = strdup(s+1);
                    }
               }
               free(info);
               return 1;
          }
          token = info_ptr;
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
     (*var)->ref = NULL;
     (*var)->alt = NULL;
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

     if (NULL == (*var)) {
          return;
     }

     free((*var)->chrom);
     free((*var)->id);
     free((*var)->ref);
     free((*var)->alt);
     free((*var)->filter);
     free((*var)->info);

     free((*var)->format);
     for (i=0; i<(*var)->num_samples; i++) {
          free((*var)->samples[i]);
     }
     free((*var)->samples);

     free(*var);
}

void vcf_cp_var(var_t **dest, var_t *src)
{
     int i;
     vcf_new_var(dest);
     (*dest)->chrom = strdup(src->chrom);
     (*dest)->pos = src->pos;
     if (src->id) {
          (*dest)->id = strdup(src->id);
     }
     if (src->ref) {
          (*dest)->ref = strdup(src->ref);
     }
     if (src->alt) {
          (*dest)->alt = strdup(src->alt);
     }
     (*dest)->qual = src->qual;
     if (src->filter) {
          (*dest)->filter = strdup(src->filter);
     }
     if (src->info) {
          (*dest)->info = strdup(src->info);
     }
     if (src->format) {
          (*dest)->format = strdup(src->format);
     }
     (*dest)->num_samples = src->num_samples;
     if (src->num_samples>0) {
          (*dest)->samples = malloc(src->num_samples * sizeof(char*));
          for (i=0; i<src->num_samples; i++) {
               (*dest)->samples[i] = strdup(src->samples[i]);
          }
     }
}

void vcf_write_var(vcf_file_t *vcf_file, const var_t *var)
{
     /* in theory all values are optional */

     vcf_printf(vcf_file, "%s\t%ld\t%s\t%s\t%s\t",
             NULL == var->chrom ? VCF_MISSING_VAL_STR : var->chrom,
             var->pos + 1,
             NULL == var->id ? VCF_MISSING_VAL_STR : var->id,
             var->ref,
             var->alt);
     if (var->qual>-1) {
          vcf_printf(vcf_file, "%d", var->qual);
     } else {
          vcf_printf(vcf_file, "%c", VCF_MISSING_VAL_CHAR);
     }

     vcf_printf(vcf_file, "\t%s\t%s",
             var->filter ? var->filter : VCF_MISSING_VAL_STR,
             var->info ? var->info : VCF_MISSING_VAL_STR);

     if (var->format) {
          int i=0;
          vcf_printf(vcf_file, "\t%s", var->format);
          for (i=0; i<var->num_samples; i++) {
               vcf_printf(vcf_file, "\t%s", var->samples[i]);
          }
     }
     vcf_printf(vcf_file, "\n");
}


char *
vcf_var_add_to_info(var_t *var, const char *info_str)
{
     if (!var || !info_str) {
          return NULL;
     }
     var->info = realloc(var->info,
                         (strlen(var->info) + strlen(info_str)
                          + 1/*;*/ + 1/*\0*/) * sizeof(char));
     if (!var->info) {
          return NULL;
     }
     if (strlen(var->info)) {
          if (0 == strcmp(var->info, VCF_MISSING_VAL_STR)) {
               var->info[0] = '\0';
          } else {
               (void) strcat(var->info, ";");
          }
     }
     (void) strcat(var->info, info_str);
     return var->info;
}

char *
vcf_var_add_to_filter(var_t *var, const char *filter_name)
{
     if (! filter_name || ! var) {
          return NULL;
     }
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
     if (! var->filter) {
          fprintf(stderr, "FATAL: couldn't allocate memory at %s:%s():%d\n",
                  __FILE__, __FUNCTION__, __LINE__);
          return NULL;
     }

     /* add */
     if (strlen(var->filter)) {
          (void) strcat(var->filter, ";");
     }
     (void) strcat(var->filter, filter_name);

     return var->filter;
}


int vcf_get_dp4(dp4_counts_t *dp4, var_t *var)
{
     const char delimiter[] = ",";
     char *token;
     char *dp4_char = NULL;
     char *dp4_char_cp;
     int i = 0;

     if ( ! vcf_var_has_info_key(&dp4_char, var, "DP4")) {
          memset(dp4, -1, sizeof(dp4_counts_t)); /* -1 = error */
          return 1;
     }
     /* note: strsep modifies ptr */
     dp4_char_cp = strdup(dp4_char);
     free(dp4_char);
     dp4_char = dp4_char_cp;

     i = 0;
     /* note: strsep modifies ptr */
     while (NULL != (token = strsep(& dp4_char, delimiter))) {
          int val = strtol(token, (char **) NULL, 10); /* = atoi */
          if (i==0) {
               dp4->ref_fw = val;
          } else if (i==1) {
               dp4->ref_rv = val;
          } else if (i==2) {
               dp4->alt_fw = val;
          } else if (i==3) {
               dp4->alt_rv = val;
          }
          i += 1;
     }
     free(dp4_char_cp);
     if (i != 4) {
          memset(dp4, -1, sizeof(dp4_counts_t)); /* -1 = error */
          return 1;
     }
     return 0;
}


/* var->info allocated here. caller has to free */
void vcf_var_sprintf_info(var_t *var,
                          const int dp, const float af, const int sb,
                          const dp4_counts_t *dp4,
                          const int indel, const int hrun, 
                          const int consvar)
{
     char buf[LINE_BUF_SIZE];
     snprintf(buf, sizeof(buf)-32, /* leave some for INDEL and other flags below */
              "DP=%d;AF=%f;SB=%d;DP4=%d,%d,%d,%d",
              dp, af, sb, dp4->ref_fw, dp4->ref_rv, dp4->alt_fw, dp4->alt_rv);
     if (indel) {
          sprintf(buf, "%s;INDEL", buf);
          sprintf(buf, "%s;HRUN=%d", buf, hrun);
     }
     if (consvar) {
          sprintf(buf, "%s;CONSVAR", buf);
     }

     var->info = strdup(buf);

     /* FIXME format and samples not supported */
}


void vcf_write_header(vcf_file_t *vcf_file, const char *header)
{
#if 0
     fprintf(stderr, "TMP DEBUG: writing header %s", header);
     fprintf(stderr, "TMP DEBUG: vcf_file path = %s\n", vcf_file->path);
     fprintf(stderr, "TMP DEBUG: vcf_file is_bgz = %d\n", vcf_file->is_bgz);
     fprintf(stderr, "TMP DEBUG: vcf_file fh = %p\n", vcf_file->fh);
     fprintf(stderr, "TMP DEBUG: vcf_file fh_bgz = %p\n", vcf_file->fh_bgz);
     fprintf(stderr, "TMP DEBUG: vcf_file mode = %c\n", vcf_file->mode);
#endif
     vcf_printf(vcf_file, "%s", header);
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

     vcf_printf(vcf_file, "##fileformat=VCFv4.0\n");
     vcf_printf(vcf_file, "##fileDate=%s\n", tbuf);
     if (src) {
          vcf_printf(vcf_file, "##source=%s\n", src);
     }
     if (reffa) {
          vcf_printf(vcf_file, "##reference=%s\n", reffa);
     }
     vcf_printf(vcf_file, "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Raw Depth\">\n");
     vcf_printf(vcf_file, "##INFO=<ID=AF,Number=1,Type=Float,Description=\"Allele Frequency\">\n");
     vcf_printf(vcf_file, "##INFO=<ID=SB,Number=1,Type=Integer,Description=\"Phred-scaled strand bias at this position\">\n");
     vcf_printf(vcf_file, "##INFO=<ID=DP4,Number=4,Type=Integer,Description=\"Counts for ref-forward bases, ref-reverse, alt-forward and alt-reverse bases\">\n");
     vcf_printf(vcf_file, "##INFO=<ID=INDEL,Number=0,Type=Flag,Description=\"Indicates that the variant is an INDEL.\">\n");
     vcf_printf(vcf_file, "##INFO=<ID=CONSVAR,Number=0,Type=Flag,Description=\"Indicates that the variant is a consensus variant (as opposed to a low frequency variant).\">\n");
     vcf_printf(vcf_file, "##INFO=<ID=HRUN,Number=1,Type=Integer,Description=\"Homopolymer length to the right of report indel position\">\n");
     vcf_printf(vcf_file, "%s\n", VCF_HEADER);
}


/* parse header, i.e. meta info until and including header from vcf
 * file. will allocate memory for header. caller has to free. returns
 * 0 on success. -1 on failure on which a minimal header is set anyway
 * and you should rewind.
 */
int vcf_parse_header(char **header, vcf_file_t *vcf_file)
{
     char line[LINE_BUF_SIZE];
     const int MAX_HEADER_LEN = 10000;
     int line_no = 0;

     /* make sure strlen below will work on header */
     (*header) = malloc(sizeof(char));
     (*header)[0] = '\0';

     line_no = 0;
     while (1) {
          char *rc = vcf_file_gets(vcf_file, sizeof(line), line);
          if (++line_no>MAX_HEADER_LEN) {
               break;
          }
          if (NULL == rc) {
               break;
          }
#if 0
          fprintf(stderr, "Got line %s\n", line);
#endif
          (*header) = realloc((*header), (strlen(*header) + strlen(line) + 1 /* '\0' */) * sizeof(char));
          (void) strcat((*header), line);
          if (strlen(line) >= strlen(VCF_HEADER)) {
               if (0 == strncmp(line, VCF_HEADER, strlen(VCF_HEADER))) {
                    return 0;
               }
          }
     }
     /* failed. set default header */
     (*header) = realloc((*header), (strlen(VCF_HEADER) + 1 + 1 /* \n+\0 */) * sizeof(char));
     (void) strcpy(*header, VCF_HEADER);
     (void) strcat(*header, "\n");

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


int vcf_parse_var_from_line(char *line, var_t *var)
{
     const char delimiter[] = "\t";
     char *token;
     char *line_ptr;
     int field_no = 0;
     char *line_backup;

     chomp(line);
     line_ptr = line;
     line_backup = strdup(line);
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
               var->ref = strdup(token);

          } else if (5 == field_no) {
               var->alt = strdup(token);

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
     if (field_no<5) {
          LOG_WARN("Parsing of variant incomplete. Only got %d fields. Need at least 5 (line=%s)\n", field_no, line_backup);
          return -1;
     }
     /* allow lenient parsing and fill in missing values*/
     if (field_no<8) {
          /* 6-8: qual, filter, info with qual already set */
          var->filter = calloc(2, sizeof(char));
          var->filter[0] = VCF_MISSING_VAL_CHAR;
          var->info = calloc(2, sizeof(char));
          var->info[0] = VCF_MISSING_VAL_CHAR;
     }

     free(line_backup);

     return 0;
}


/* parse one variant from stream. returns -1 on error or EOF.
 * note, multi-allelic entries are not treated specially. returns non-null on error
 */
int vcf_parse_var(vcf_file_t *vcf_file, var_t *var)
{
     char line[LINE_BUF_SIZE];
     char *rc;

     rc = vcf_file_gets(vcf_file, sizeof(line), line);
     if (NULL == rc) {
          return -1;
     }
     return vcf_parse_var_from_line(line, var);
}


/* parse all variants from stream and return number of parsed vars or
 * -1 on error. memory for vars will be allocated here.
 */
int vcf_parse_vars(var_t ***vars, vcf_file_t *vcf_file, int only_passed)
{
     int rc;
     int num_vars = 0;

     (*vars) = malloc(1 * sizeof(var_t*));

     while (1) { 
          var_t *var;
          vcf_new_var(&var);
          rc = vcf_parse_var(vcf_file, var);
          if (rc) {
               /* would be nice to distinguish between eof and error */
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

     /* make sure to insert before VCF_HEADER */

     token = strstr(*header, VCF_HEADER);
     if (! token) {
          LOG_WARN("%s\n", "Can't add info to empty header, because header line is missing");
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
     (void) strcat(*header, VCF_HEADER);
     (void) strcat(*header, "\n");
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
     LOG_INFO("Using %s (%s gzipped)\n", path_in, gzip_in ? "is" : "not");
     if (vcf_file_open(& vcf_file_in, path_in, gzip_in, 'r')) {
          LOG_FATAL("%s\n", "vcf_file_open() failed");
          exit(1);
     }

     if (HAS_GZIP_EXT(path_out)) {
          gzip_out = 1;
     } else {
          gzip_out = 0;
     }
     LOG_INFO("Using %s (%s gzipped)\n", path_out, gzip_out ? "is" : "not");
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


     num_vars = vcf_parse_vars(& vcf_file_in, &vars);
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
