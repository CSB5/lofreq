#include <ctype.h>
#include <stdio.h>
#include <time.h>

#include "faidx.h"
#include "sam.h" 
#include "log.h"
#include "lofreq_indel_quality.h"

char DINDELQ[] = "!MMMLKEC@=<;:988776"; // 1-based 18
char DINDELQ2[] = "!CCCBA;963210/----,"; // *10

typedef struct {
     samfile_t *in;
     bamFile out;
     uint8_t defq[1000];
} tmpstruct_t_default;

typedef struct {
     samfile_t *in;
     bamFile out;
     faidx_t *fai;
     int *hpcount;
     int rlen;
     uint32_t tid;
} tmpstruct_t_dindel;

#define prob_to_sangerq(p) (p < 0.0 + DBL_EPSILON ? 126+1 : ((int)(-10 * log10(p))+33))
#define encode_q(q) (uint8_t)(q < 33 ? '!' : (q > 126 ? '~' : q))

static int default_fetch_func(bam1_t *b, void *data)
{
     tmpstruct_t_default *tmp = (tmpstruct_t_default*)data;
     bam1_core_t *c = &b->core;

     uint8_t indelq[c->l_qseq+1];

     memcpy(&indelq, &tmp->defq, c->l_qseq);
     indelq[c->l_qseq] = '\0';

     uint8_t *to_delete;
     to_delete = bam_aux_get(b, "BI");
     if (to_delete) bam_aux_del(b, to_delete);
     bam_aux_append(b, "BI", 'Z', c->l_qseq+1, indelq);


     to_delete = bam_aux_get(b, "BD");
     if (to_delete) bam_aux_del(b, to_delete);
     bam_aux_append(b, "BD", 'Z', c->l_qseq+1, indelq);

     bam_write1(tmp->out, b);
     return 0;
}

/* Stores an array of ints that corresponds to the length of the homopolymer 
 * at the start of each homopolymer*/
int find_homopolymers(char *query, int *count, int qlen)
{
     int i, j;
     int curr_i = 0;
     int curr_count = 1;
     for (i = 1; i < qlen; i++) {
          if (query[i] == query[curr_i]) {
               curr_count += 1; // keep incrementing count if in homopolymer region
          } else {
               count[curr_i] = curr_count; // record length of homopolymer region
               for (j = curr_i+1; j < i; j++) {
                    count[j] = 1; // all other positions get a count of 1
               }
               curr_i = i;
               curr_count = 1;
          }
     }
     if (curr_i < i) { // take care of edge case at the end of the read
          count[curr_i] = curr_count;
          for (j = curr_i+1; j < i; j++) {
               count[j] = 1;
          }
     }
     return 0;
}

static int dindel_fetch_func(bam1_t *b, void *data)
{
     tmpstruct_t_dindel *tmp = (tmpstruct_t_dindel*)data;
     bam1_core_t *c = &b->core;
     int rlen;

     // skip all reads that are not properly paired
     if (! (c->flag & BAM_FPROPER_PAIR)) {
          // fprintf(stderr, "skipping read: %s at pos %d\n", bam1_qname(b), c->pos);
          return 0;
     }

     // get the reference sequence and compute homopolymer array
     if (tmp->tid != c->tid) {
          //fprintf(stderr, "fetching reference sequence %s\n",
          //                tmp->in->header->target_name[c->tid]);
          char *ref = fai_fetch(tmp->fai, tmp->in->header->target_name[c->tid], &rlen);
          int rlen = strlen(ref);
          tmp->tid = c->tid;
          if (tmp->hpcount) free(tmp->hpcount);
          tmp->hpcount = (int*)malloc(rlen*sizeof(int));
          find_homopolymers(ref, tmp->hpcount, rlen);
          free(ref);
          tmp->rlen = rlen;
          //fprintf(stderr, "fetched reference sequence\n");
     }

     // parse the cigar string
     uint32_t *cigar = bam1_cigar(b);
     uint8_t indelq[c->l_qseq+1];
     //fprintf(stderr, "l_qseq:%d\n", c->l_qseq);
     int i;
     int x = c->pos; // coordinate on reference
     int y = 0; // coordinate on query
     for (i = 0; i < c->n_cigar; ++i) {
          int j, oplen = cigar[i]>>4, op = cigar[i]&0xf;
          if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF) {
               for (j = 0; j < oplen; j++) {
                    //fprintf(stderr, "query:%d, ref:%d, count:%d\n", 
                    //        y, x, tmp->hpcount[x+1]);
                    indelq[y] = (x > tmp->rlen-2) ? DINDELQ[0] : (tmp->hpcount[x+1]>18? 
                         DINDELQ[0] : DINDELQ[tmp->hpcount[x+1]]);
                    x++; 
                    y++;
               }
          } else if (op == BAM_CHARD_CLIP) { // do nothing
          } else if (op == BAM_CDEL) {
               x += oplen;
          } else if (op == BAM_CINS || op == BAM_CSOFT_CLIP) { 
               for (j = 0; j < oplen; j++) {
                    //fprintf(stderr, "query:%d, ref:%d\n", y, x);
                    indelq[y] = DINDELQ[0];
                    y++;
               }
          } else {
               exit(1);
          }
     }
     indelq[y] = '\0';

     bam_aux_append(b, "BI", 'Z', c->l_qseq+1, indelq);
     bam_aux_append(b, "BD", 'Z', c->l_qseq+1, indelq);

     bam_write1(tmp->out, b);
     return 0;
}

int add_default(char *argv[])
{
	tmpstruct_t_default tmp;
	if ((tmp.in = samopen(argv[3], "rb", 0)) == 0) {
              fprintf(stderr, "default_indel_quality: Failed to open BAM file %s\n", argv[3]);
              return 1;
         }

         uint8_t q = encode_q(atoi(argv[4])+33);

         int i;
         for (i = 0; i < 999; i++) {
              tmp.defq[i] = q;
         }
         tmp.defq[i] = '\0';

         tmp.out = bam_dopen(fileno(stdout), "w");
         bam_header_write(tmp.out, tmp.in->header);
          
         bam1_t *b = bam_init1();
         int count = 0;
         while (samread(tmp.in, b) >= 0) {
              count++;
              default_fetch_func(b, &tmp); 
         }
         bam_destroy1(b);
        
         samclose(tmp.in);
         bam_close(tmp.out);
	 fprintf(stderr, "Processed %d reads\n", count);
	 return 0;
}

int add_dindel(char *argv[]){
	tmpstruct_t_dindel tmp;
	
	if ((tmp.in = samopen(argv[3], "rb", 0)) == 0) {
             fprintf(stderr, "indel_quality: Failed to open BAM file %s\n", argv[1]);
             return 1;
        }
        if ((tmp.fai = fai_load(argv[4])) == 0) {
             fprintf(stderr, "indel_quality: Failed to open .fa file%s\n", argv[2]);
             return 1;
        }

        tmp.out = bam_dopen(fileno(stdout), "w");
        bam_header_write(tmp.out, tmp.in->header);
        
        bam1_t *b = bam_init1();
        tmp.tid = -1;
        tmp.hpcount = 0;
        tmp.rlen = 0;
        int count = 0;
        fprintf(stderr, "WARNING: DO NOT REALIGN AFTER CALIBRATION\n");
        while (samread(tmp.in, b) >= 0) {
             count++;
             dindel_fetch_func(b, &tmp); 
        }
        bam_destroy1(b);
        
        if (tmp.hpcount) free(tmp.hpcount);
        samclose(tmp.in);
        bam_close(tmp.out);
        fai_destroy(tmp.fai);
	fprintf(stderr, "Processed %d reads\n", count);
	return 0;
}

int main_indel_quality(int argc, char *argv[])
{
     
     if (argc < 5) {
          fprintf(stderr, "Usage: lofreq <main_flag> <subflag> <in.bam> <defaultq/ref.fa>\n");
          return 1;
     }
     if (argc == 5) {
          
          time_t now = time(NULL);
          char date[100];
          strftime(date, 100, "%c", localtime(&now));
          fprintf(stderr, "indel_quality: Started at %s.\n", date);
	
	  if (strcmp(argv[2], "default")==0){
		// adding default values
		return add_default(argv);
	  } else if (strcmp(argv[2], "dindel") == 0){
		//adding dindel values
		return add_dindel(argv);
	  } else{
		LOG_FATAL("Unrecognized command '%s'\n", argv[2]);
	  	return 1;
	  }

	  now = time(NULL);
          strftime(date, 100, "%c", localtime(&now));
          fprintf(stderr, "Ended at %s.\n", date);

     } else {
          fprintf(stderr, "Usage: lofreq <main_flag> <subflag> <in.bam> <defaultq/ref.fa>\n");
          return 1;
     }
     return 0;
}
