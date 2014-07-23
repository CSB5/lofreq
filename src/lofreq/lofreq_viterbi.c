#include <ctype.h>
#include <stdio.h>
#include <time.h>
#include <getopt.h>
#include <stdlib.h>

#include "faidx.h"
#include "sam.h"
#include "viterbi.h"
#include "log.h"
/* added lofreq_viterbi.h" */
#include "lofreq_viterbi.h"
#include "utils.h"

#define SANGERQUAL_TO_PHRED(c) ((int)(c)-33)

// FIXME:
//	viterbi output not necessarily sorted anymore!!
//	remedy: pipe through samtools sort - (warning added to terminal output)
//
//	auto clipping of Q2 tails (complicated stuff!)

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


// function checks if alignment is made of all Q2s
// if not, returns remaining values so that median 
// can be calculated
int check_Q2(char *bqual, int *num){
	//printf("inside check\n");
	int is_all_Q2 = 1;
	int i, pom = 0;
	int l = strlen(bqual);
//	char c;
	*num = 0;
	//printf("Length %d\n", l);
	for (i=0; i<l; i++){
		if (SANGERQUAL_TO_PHRED(bqual[i]) != 2){
			pom++;
			is_all_Q2 = 0;
			//printf("Found one without Q2, num is %d\n", pom);
		//	while(!(c=fgetc(stdin)));
		}
	}		
	*num = pom;
	return is_all_Q2;
}

void remain(char *bqual, int *remaining){
	int pom = 0;
	int i, q;
//	char c;
	int l = strlen(bqual);
	for (i=0; i<l; i++){
		q = SANGERQUAL_TO_PHRED(bqual[i]);
		if (q != 2){
			remaining[pom] = q;
		//	printf("Remaining[pom] %d\n",remaining[pom]);
		//	while(!(c=fgetc(stdin)));
			pom++;
		}
	}	
}

static int fetch_func(bam1_t *b, void *data, int flag, int quality)
{
     // see https://github.com/lh3/bwa/blob/426e54740ca2b9b08e013f28560d01a570a0ab15/ksw.c 
     // for optimizations and speedups

	//char chr;
     tmpstruct_t *tmp = (tmpstruct_t*)data;
     bam1_core_t *c = &b->core;
     uint8_t *seq = bam1_seq(b);
     uint32_t *cigar = bam1_cigar(b);
     int reflen;
	
	 // FIXME: removal of unupdated flags
     // IDEA (Andreas) : flags (NM, MC, MD, AS) might become invalid
     // removal is optional (user flag)
     if (flag == 1){
	uint8_t *old_nm = bam_aux_get(b, "NM");
	if (old_nm) bam_aux_del(b, old_nm);

	uint8_t *old_mc = bam_aux_get(b, "MC");
	if (old_mc) bam_aux_del(b, old_mc);

	uint8_t *old_md = bam_aux_get(b, "MD");
	if (old_md) bam_aux_del(b, old_md);

	uint8_t *old_as = bam_aux_get(b, "AS");
	if (old_as) bam_aux_del(b, old_as);
     }

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
	FILE *f;
	f = fopen("dogadaji.txt","w"); 
	int len_remaining = 0;
	if (check_Q2(bqual, &len_remaining)){
	//	printf("All are Q2\n");
	//	while(!(chr=fgetc(stdin)));
		bam_write1(tmp->out, b);
		fprintf(f,"All Q2\n");
        return 0;
     }
	fprintf(f,"Continuing\n");
//	printf("starting\n");
//	printf("len_remaining %d\n",len_remaining);
//	while(!(chr=fgetc(stdin)));
	int remaining[len_remaining+1];
	remain(bqual, remaining);
	remaining[len_remaining] = '\0';
	if (quality < 0){
		quality = int_median(remaining, len_remaining);
		fprintf(f, "New quality %d\n", quality);
//		printf("New quality %d\n", quality);
//		while(!(chr=fgetc(stdin)));
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
     int shift = viterbi(ref, query, bqual, aln, quality);
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
	 fclose(f);
     return 0;
}

/* renamed the main function */
int main_viterbi(int argc, char *argv[])
{
     tmpstruct_t tmp;
     static int usrflg = 0;
     static int quality = -1;
     int c;
    
	time_t now = time(NULL);
	char date[100];
	strftime(date, 100, "%c", localtime(&now));
	fprintf(stderr, "started at %s.\n", date);

     /*should be 2 args, lofreq and viterbi */
     if (argc == 2) {
          fprintf(stderr, "Usage: lofreq viterbi [options] in.bam\n");
	  fprintf(stderr, "Options:\n");
	  fprintf(stderr, "	Reference\n");
	  fprintf(stderr, "	    -f | --ref FILE	Indexed reference fasta file [null]\n");
	  fprintf(stderr, "	Flags\n");
	  fprintf(stderr, "	    -d | --delete	Deletes flags (MC, MD, NM, AS) from BAM file\n");
	  fprintf(stderr, "	Q2 qualities\n");
	  fprintf(stderr, "	    -q | --qualities Q	Replaces Q2 qualities with given quality Q\n");

	  fprintf(stderr, "WARNING: 	BAM file is unsorted after viterbi, pipelining through\n");
	  fprintf(stderr, "		'samtools sort -' obligatory\n");
          return 1;
     }

		//make structure for long options
		static struct option long_options[] = {
			{"ref", required_argument, NULL, 'f'},
			{"delete", no_argument, NULL, 'd'},
			{"qualities", required_argument, NULL, 'q'},
			{0,0,0,0}
		};
		
		// make string of arguments
		static const char *long_opts_str = "df:q:";
		
		//getopt_long stores option index
		int long_option_index = 0;
		
		while((c = getopt_long(argc-1, argv+1, long_opts_str, long_options, &long_option_index)) != -1){
			switch (c){
				case 'f':
					if(access(optarg,F_OK)== -1){
						LOG_FATAL("Reference fasta file %s does not exist. Exiting...\n", optarg);
						return 1;
				}
				tmp.fai = fai_load(optarg);	
				break;
			case 'd':
				usrflg = 1;
				break;
			case 'q':
				quality = atoi(optarg);
				fprintf(stderr, "Given quality %d \n", atoi(optarg));
				break;
			case '?':
				LOG_FATAL("%s\n", "Unrecognized arguments found. Exiting\n");
				break;
			default:
				break;
		}
		}
		
		// get bam file
		if (1 != argc-optind-1){
			fprintf(stderr, "Need exactly one BAM file as last argument\n");
			return 1;
		}
		if ((tmp.in = samopen((argv+optind+1)[0], "rb",0)) == 0){
			fprintf(stderr, "Failed to open BAM file %s. Exiting...\n", (argv+optind+1)[0]);
			return 1;
		}

		tmp.out = bam_dopen(fileno(stdout), "w");
		bam_header_write(tmp.out, tmp.in->header);

		bam1_t *b = bam_init1();
		tmp.tid = -1;
		tmp.ref = 0;
		while (samread(tmp.in, b) >= 0){
			fetch_func(b, &tmp, usrflg, quality);
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
}
