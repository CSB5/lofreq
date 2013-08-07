#!/bin/bash

# Test that we get the number of expected SNVs on the pseudo-clonal data-set

source lib.sh || exit 1


basedir=data/ecoli-clone/
bam=$basedir/Ecoli_K12_MG1655_NC_000913_bwa-sampe-unique-rg_recal.bam
reffa=$basedir/Ecoli_K12_MG1655_NC_000913.fa.gz
#truesnv=$basedir/denv2-pseudoclonal_true-snp.vcf

outdir=$(mktemp -d -t $(basename $0).XXXXXX)
outvcf=$outdir/$(basename $bam .bam).vcf
log=$outdir/log.txt

KEEP_TMP=0

cmd="$LOFREQ call -f $reffa -E -o $outvcf $bam"
echodebug "cmd=$cmd"
if ! eval $cmd >> $log 2>&1; then
    echoerror "The following command failed (see $log for more): $cmd"
    exit 1
fi


MAX_SNVS=20
nsnvs=$(grep -c '^[^#]' $outvcf)
if [ $nsnvs -ge $MAX_SNVS ]; then
    echoerror "Expected less then $MAX_SNVS on this clonal dataset"
    exit 1
else
    echook "Got $nsnvs SNVs for this clonal dataset which is okay (below limit of $MAX_SNVS)"
fi


if [ $KEEP_TMP -eq 1 ]; then
    echowarn "Not deleting tmp dir $outdir"
else 
    rm $outdir/*
    rmdir $outdir
fi
