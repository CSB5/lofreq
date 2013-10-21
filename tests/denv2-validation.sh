#!/bin/bash

# Test that we get the number of expected SNVs on the pseudo-clonal data-set

source lib.sh || exit 1

basedir=./data/denv2-dpcr-validated
bam1=$basedir/CTTGTA_2_remap_razers-i92_peakrem_corr.bam
bam2=$basedir/GGCTAC_2_remap_razers-i92_peakrem_corr.bam
reffa=$basedir/consensus.fa

outdir=$(mktemp -d -t $(basename $0).XXXXXX)
vcfout1=$outdir/$(basename $bam1 .bam).vcf
vcfout2=$outdir/$(basename $bam2 .bam).vcf
vcfinter=$outdir/intersection.vcf

log=$outdir/log.txt

KEEP_TMP=0

cmd="$LOFREQ call -B -f $reffa -o $vcfout1 $bam1"
if ! eval $cmd >> $log 2>&1; then
    echoerror "The following command failed (see $log for more): $cmd"
    exit 1
fi

cmd="$LOFREQ call -B -f $reffa -o $vcfout2 $bam2"
if ! eval $cmd >> $log 2>&1; then
    echoerror "The following command failed (see $log for more): $cmd"
    exit 1
fi

cmd="$LOFREQ vcfset -a intersect -1 $vcfout1 -2 $vcfout2 -o $vcfinter"
if ! eval $cmd >> $log 2>&1; then
    echoerror "The following command failed (see $log for more): $cmd"
    exit 1
fi

N_PRESENT=7
n_present=$(for pos in 5914 6843 598 5025 1687 9941 4828; do grep "^consensus[^0-9]*$pos" $vcfinter; done | wc -l)
N_ABSENT=0
n_absent=$(for pos in 7035 7404; do grep "^consensus[^0-9]*$pos" $vcfinter; done | wc -l)

if [ $n_present -ne $N_PRESENT ]; then
    echoerror "Expected $N_PRESENT but got $n_present SNVs"
    exit 1
fi

if [ $n_absent -ne $N_ABSENT ]; then
    echoerror "Expected $N_ABSENT but got $n_absent SNVs"
    exit 1
fi

echook "Got expected number of present/absent SNVs"

if [ $KEEP_TMP -eq 1 ]; then
    echowarn "Not deleting tmp dir $outdir"
else 
    rm  $outdir/*
    rmdir $outdir
fi

