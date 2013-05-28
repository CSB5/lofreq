#!/bin/bash

# Test that we get the number of expected SNVs on the pseudo-clonal data-set

source lib.sh || exit 1


basedir=data/denv2-pseudoclonal
bam=$basedir/denv2-pseudoclonal.bam
reffa=$basedir/denv2-pseudoclonal_cons.fa
bed=$basedir/denv2-pseudoclonal_incl.bed
#truesnv=$basedir/denv2-pseudoclonal_true-snp.vcf

outdir=$(mktemp -d -t $(basename $0).XXXXXX)
outraw_nobaq=$outdir/raw_nobaq.vcf
outraw_baq=$outdir/raw_baq.vcf
log=$outdir/log.txt

KEEP_TMP=0

cmd="$LOFREQ call -f $reffa -l $bed -o $outraw_nobaq $bam"
if ! eval $cmd >> $log 2>&1; then
    echoerror "The following command failed (see $log for more): $cmd"
    exit 1
fi

cmd="$LOFREQ call -E -f $reffa -l $bed -o $outraw_baq $bam"
if ! eval $cmd >> $log 2>&1; then
    echoerror "The following command failed (see $log for more): $cmd"
    exit 1
fi


ndiff=$($LOFREQ vcfset -a complement -1 $outraw_nobaq -2 $outraw_baq  | grep -c '^[^#]')
if [ $ndiff -lt 1 ]; then
    echoerror "Expected more SNVs with BAQ switched off (check $outraw_nobaq and $outraw_baq)"
    exit 1
else
    echook "Got $ndiff more SNVs if BAQ is off"
fi



if [ $KEEP_TMP -eq 1 ]; then
    echowarn "Not deleting tmp dir $outdir"
else 
    rm  $outdir/*
    rmdir $outdir
fi

