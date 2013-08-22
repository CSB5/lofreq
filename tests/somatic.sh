#!/bin/bash

# FIXME:add-doc

source lib.sh || exit 1


KEEP=1

bam_n=./data/somatic/CHH966-normal-100x-100pur-hg19.chr22-bed-only.bam
bam_t=./data/somatic/CHH966-tumor-100x-10pur-hg19.chr22-bed-only.bam
bed=./data/somatic/SeqCap_EZ_Exome_v3_primary_lib_extend_no_overlap_minus300.chr22.bed
reffa=./data/somatic/hg19_chr22.fa.gz
truesnv=./data/somatic/hg19_chr22_true_snv.vcf
#outprefix=$(mktemp -t $(basename $0).XXXXXX)
outprefix=$(mktemp -t $(basename $0));#XXXXXX needed on linux?
finalout=${outprefix}lofreq_somatic_final.vcf
#$LOFREQ somatic -n $bam_n -t $bam_t -f $reffa -l $bed -o $outprefix || exit 1
#$LOFREQ somatic -S 10 -n $bam_n -t $bam_t -f $reffa -l $bed -o $outprefix || exit 1
$LOFREQ somatic -F 0.01 -n $bam_n -t $bam_t -f $reffa -l $bed -o $outprefix || exit 1
n_intersect=$($LOFREQ vcfset -1 $truesnv -2 $finalout -a intersect | grep -vc '^#')
if [ "$n_intersect" -lt 2 ]; then
	echoerror "Expected at least two true predictions but got $n_intersect (compare $finalout and $truesnv)"
else
	echook "Got $n_intersect true predictions"
    if [ $KEEP -eq 1 ]; then
        echodebug "Not deleting ${outprefix}*vcf*"
    else
	    rm ${outprefix}*vcf*
    fi
fi



