/* -*- c-file-style: "k&r" -*-
 *
 *
 * FIXME missing license
 *
 */
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "vcf.h"

#define MISSING_VAL_STR "."
#define MISSING_VAL_CHAR '.'


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
}


void vcf_free_var(var_t **var)
{
     free((*var)->chrom);
     free((*var)->id);
     free((*var)->filter);
     free((*var)->info);
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
}


/* var->info allocated here. caller has to free */
void vcf_var_sprintf_info(var_t *var,
                         const int *dp, const float *af, const int *sb,
                         const dp4_counts_t *dp4,
                         const int indel, const int consvar)
{
     char buf[1024];          
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
}


void vcf_write_header(FILE *stream, const char *srcprog, const char *reffa)
{
     char tbuf[9];
     struct tm tm;
     time_t t;

     t = time(0);
     localtime_r(&t, &tm);
     strftime(tbuf, 9, "%Y%m%d", &tm);
     
     fprintf(stream, "##fileformat=VCFv4.1\n");
     fprintf(stream, "##fileDate=%s\n", tbuf);
     fprintf(stream, "##source=\"%s\"\n", srcprog);
     fprintf(stream, "##reference=%s\n", reffa);

     fprintf(stream, "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Raw Depth\">\n");
     fprintf(stream, "##INFO=<ID=AF,Number=1,Type=Float,Description=\"Allele Frequency\">\n");    
     fprintf(stream, "##INFO=<ID=SB,Number=1,Type=Integer,Description=\"Phred-scaled strand bias at this position\">\n");    

     fprintf(stream, "##INFO=<ID=DP4,Number=4,Type=Integer,Description=\"Counts for ref-forward bases, ref-reverse, alt-forward and alt-reverse bases\">\n");

     fprintf(stream, "##INFO=<ID=INDEL,Number=0,Type=Flag,Description=\"Indicates that the variant is an INDEL.\">\n");
     fprintf(stream, "##INFO=<ID=CONSVAR,Number=0,Type=Flag,Description=\"Indicates that the variant is a consensus variant (as opposed to a low frequency variant).\">\n");

     fprintf(stream, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");
}
