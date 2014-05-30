#!/bin/bash

# Test that we get the number of expected SNVs on the pseudo-clonal data-set

source lib.sh || exit 1


basedir=data/denv2-pseudoclonal
bam=$basedir/denv2-pseudoclonal.bam
reffa=$basedir/denv2-pseudoclonal_cons.fa
bed=$basedir/denv2-pseudoclonal_incl.bed
truesnv=$basedir/denv2-pseudoclonal_true-snp.vcf

outdir=$(mktemp -d -t $(basename $0).XXXXXX)
outraw_def=$outdir/raw_def.vcf
outfinal_def=$outdir/final_def.vcf
log=$outdir/log.txt

KEEP_TMP=0

cmd="$LOFREQ call -b dynamic -f $reffa -l $bed -o $outraw_def $bam"
if ! eval $cmd >> $log 2>&1; then
    echoerror "The following command failed (see $log for more): $cmd"
    exit 1
fi
cmd="$LOFREQ filter --only-passed -i $outraw_def -o $outfinal_def"
if ! eval $cmd >> $log 2>&1; then
    echoerror "The following command failed (see $log for more): $cmd"
    exit 1
fi


ndiff=$($LOFREQ vcfset -a complement --only-passed -1 $outfinal_def -2 $truesnv  | grep -c '^[^#]')
if [ $ndiff -ne 0 ]; then
    echoerror "Found FP SNVs (not part of the list of true SNVs). Check $outdir"
    exit 1
fi

ndiff=$($LOFREQ vcfset -a intersect --only-passed -1 $outfinal_def -2 $truesnv  | grep -c '^[^#]')
nexp=229
if [ $ndiff -lt $nexp ]; then
    echoerror "Expected $nexp TP SNVs but got $ndiff. Check $outdir"
    exit 1
fi


echook "Tests passed"

if [ $KEEP_TMP -eq 1 ]; then
    echowarn "Not deleting tmp dir $outdir"
else 
    rm  $outdir/*
    rmdir $outdir
fi

