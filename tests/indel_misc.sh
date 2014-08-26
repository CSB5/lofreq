#!/bin/bash

source lib.sh || exit 1

KEEP_TMP=0
REF=data/denv2-dpcr-validated/consensus.fa
outdir=$(mktemp -d -t $(basename $0).XXXXXX)
BAM=data/denv2-dpcr-validated/CTTGTA_2_remap_razers-i92_peakrem_corr.bam 


log=$outdir/log.txt
vcf=$outdir/out.vcf
cmd="$LOFREQ call --no-default-filter --only-indels -f $REF -o $vcf $BAM"
if ! eval $cmd > $log 2>&1; then
    echoerror "LoFreq failed. Check logfile $log. Command was $cmd"
    exit 1
fi

num_indels=$(grep -vc '^#' $vcf)
if [ $num_indels -ne 0 ]; then
    echoerror "Got indels in indel free bam. See $vcf"[B
    exit 1
else
    echook "Got no indels from indel free bam."
fi

if [ $KEEP_TMP -ne 1 ]; then
	test -d $outdir && rm -rf $outdir
fi

# FIXME call on ecoli as well or test there