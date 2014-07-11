#include <ctype.h>
#include <stdio.h>
#include <time.h>

#include "faidx.h"
#include "sam.h"
#include "viterbi.h"

#define RWIN 10

typedef struct {
     samfile_t *in;
     bamFile out;
     faidx_t *fai;
     uint32_t tid;
     char *ref;
     int reflen;
} tmpstruct_t;

static void replace_cigar(bam1_t *b, int n, uint32_t *cigar)
{
	if (n != b->core.n_cigar) {
		int o = b->core.l_qname + b->core.n_cigar * 4;
		if (b->data_len + (n - b->core.n_cigar) * 4 > b->m_data) {
			b->m_data = b->data_len + (n - b->core.n_cigar) * 4;
			kroundup32(b->m_data);
			b->data = (uint8_t*)realloc(b->data, b->m_data);
		}
		memmove(b->data + b->core.l_qname + n * 4, b->data + o, b->data_len - o);
		memcpy(b->data + b->core.l_qname, cigar, n * 4);
		b->data_len += (n - b->core.n_cigar) * 4;
		b->core.n_cigar = n;
	} else memcpy(b->data + b->core.l_qname, cigar, n * 4);
}

static int fetch_func(bam1_t *b, void *data)
{
     // see https://github.com/lh3/bwa/blob/426e54740ca2b9b08e013f28560d01a570a0ab15/ksw.c 
     // for optimizations and speedups

     tmpstruct_t *tmp = (tmpstruct_t*)data;
     bam1_core_t *c = &b->core;
     uint8_t *seq = bam1_seq(b);
     uint32_t *cigar = bam1_cigar(b);
     int reflen;

     if (c->flag & BAM_FUNMAP) {
          bam_write1(tmp->out, b);
          return 0;
     }

     // fetch reference sequence if incorrect tid
     if (tmp->tid != c->tid) {
          if (tmp->ref) free(tmp->ref);
          if ((tmp->ref = 
               fai_fetch(tmp->fai, tmp->in->header->target_name[c->tid], &reflen)) == 0) {
               fprintf(stderr, "failed to find reference sequence %s\n", 
                                tmp->in->header->target_name[c->tid]);
          }
          tmp->tid = c->tid;
          tmp->reflen = reflen;
     }
     int i;

     // remove soft clipped bases
     char query[c->l_qseq+1];
     char bqual[c->l_qseq+1];

     int x = c->pos; // coordinate on reference
     int y = 0; // coordinate on query
     int z = 0; // coordinate on query w/o softclip

     int indels = 0;

     // parse cigar string
     for (i = 0; i < c->n_cigar; ++i) {
          int j, oplen = cigar[i] >> 4, op = cigar[i]&0xf;
          if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF) {
               for (j = 0; j < oplen; j++) {
                    query[z] = bam_nt16_rev_table[bam1_seqi(seq, y)];
                    bqual[z] = (char)bam1_qual(b)[y]+33;
                    x++;
                    y++;
                    z++;
               }
          } else if (op == BAM_CHARD_CLIP) {
               return 1;
          } else if (op == BAM_CDEL) {
               x += oplen;
               indels += 1;
          } else if (op == BAM_CINS) {
               for (j = 0; j < oplen; j++) {
                    query[z] = bam_nt16_rev_table[bam1_seqi(seq, y)];
                    bqual[z] = (char)bam1_qual(b)[y]+33;
                    y++;
                    z++;
               }
               indels += 1;
          } else if (op == BAM_CSOFT_CLIP) {
               for (j = 0; j < oplen; j++) {
                    y++;
               }
          } else {
               return 1;
          }
     }
     query[z] = bqual[z] = '\0';

     if (indels == 0) {
          bam_write1(tmp->out, b);
          return 0;
     }
    
     // get reference with RWIN padding
     char ref[c->l_qseq+1+indels+RWIN*2];
     int lower = c->pos - RWIN;
     lower = lower < 0? 0: lower;
     int upper = x + RWIN;
     upper = upper > tmp->reflen? tmp->reflen: upper;
     for (z = 0, i = lower; i < upper; z++, i++) {
          ref[z] = toupper(tmp->ref[i]);
     }
     ref[z] = '\0';

     // run viterbi
     char *aln = malloc(sizeof(char)*(2*(c->l_qseq)));
     int shift = viterbi(ref, query, bqual, aln);
     // fprintf(stderr, "Realigned read %s with shift %d\n", 
     //                 bam1_qname(b), shift-(c->pos-lower));

     // convert to cigar
     uint32_t *realn_cigar = 0;
     int realn_n_cigar = 0;
     
     // check if soft-clipped in the front
     int curr_oplen = cigar[0] >> 4; 
     int curr_op = cigar[0]&0xf;
     if (curr_op == BAM_CSOFT_CLIP) {
          realn_cigar = realloc(realn_cigar, (realn_n_cigar+1)*sizeof(uint32_t));
          realn_cigar[realn_n_cigar] = curr_oplen<<4 | curr_op;
          realn_n_cigar += 1;
     }
     
     // get the cigar of the realigned query
     curr_op = aln[0] == 'M' ? 0 : (aln[0] == 'I'? 1 : 2);
     curr_oplen = 1;
     for (i = 1; i < strlen(aln); i++) {
          int this_op = aln[i] == 'M' ? 0 : (aln[i] == 'I' ? 1 : 2);
          if (this_op != curr_op) {
               realn_cigar = realloc(realn_cigar, (realn_n_cigar+1)*sizeof(uint32_t));
               realn_cigar[realn_n_cigar] = curr_oplen<<4 | curr_op;
               realn_n_cigar += 1;
               curr_op = this_op;
               curr_oplen = 1;
          } else {
               curr_oplen += 1;
          }
     }
     realn_cigar = realloc(realn_cigar, (realn_n_cigar+1)*sizeof(uint32_t));
     realn_cigar[realn_n_cigar] = curr_oplen<<4 | curr_op;
     realn_n_cigar += 1; 
    
     // check if soft-clipped in the back
     curr_oplen = cigar[c->n_cigar-1] >> 4; 
     curr_op = cigar[c->n_cigar-1]&0xf;
     if (curr_op == BAM_CSOFT_CLIP) {
          realn_cigar = realloc(realn_cigar, (realn_n_cigar+1)*sizeof(uint32_t));
          realn_cigar[realn_n_cigar] = curr_oplen<<4 | curr_op;
          realn_n_cigar += 1;
     }

#if 0
     int j;
     for (j = 0; j < realn_n_cigar; j++) {
          curr_oplen = realn_cigar[j] >> 4;
          curr_op = realn_cigar[j]&0xf;
          fprintf(stderr, "op:%d oplen:%d\n", curr_op, curr_oplen);
     }
#endif

     // check if read was shifted
     if (shift-(c->pos-lower) != 0) {
          fprintf(stderr, "Read %s with shift of %d at original pos %s:%d\n", 
                           bam1_qname(b), shift-(c->pos-lower),
                           tmp->in->header->target_name[c->tid], c->pos);
          c->pos = c->pos + (shift - (c->pos - lower));
     }
     
     replace_cigar(b, realn_n_cigar, realn_cigar);
     bam_write1(tmp->out, b);

     free(aln);
     free(realn_cigar);
     return 0;
}

int main(int argc, char *argv[])
{
     tmpstruct_t tmp;

     if (argc == 1) {
          fprintf(stderr, "Usage: viterbi_realigner <in.bam> <ref.fa>\n");
          return 1;
     }
     if (argc == 3) {
         
          time_t now = time(NULL);
          char date[100];
          strftime(date, 100, "%c", localtime(&now));
          fprintf(stderr, "Started at %s.\n", date);

          if ((tmp.in = samopen(argv[1], "rb", 0)) == 0) {
               fprintf(stderr, "viterbi_realigner: Failed to open BAM file %s\n", argv[1]);
               return 1;
          }
           
          if ((tmp.fai = fai_load(argv[2])) == 0) {
               fprintf(stderr, "viterbi_realigner: Failed to open .fa file %s\n", argv[2]);
               return 1;
          }

          tmp.out = bam_dopen(fileno(stdout), "w");
          bam_header_write(tmp.out, tmp.in->header);

          bam1_t *b = bam_init1(); 
          tmp.tid = -1;
          tmp.ref = 0;
          while (samread(tmp.in, b) >= 0) {
               fetch_func(b, &tmp);
          }
          bam_destroy1(b);

          fprintf(stderr, "Done reading all reads\n");

          samclose(tmp.in);
          bam_close(tmp.out);
          if (tmp.ref) free(tmp.ref);
          fai_destroy(tmp.fai);

          now = time(NULL);
          strftime(date, 100, "%c", localtime(&now));
          fprintf(stderr, "Ended at %s.\n", date);

          return 0;
     
     } else {
          fprintf(stderr, "Usage: viterbi_realigner <in.bam> <ref.fa>\n");
          return 1;
     }

}
