#!/bin/bash

# Test that we get the number of expected SNVs on the pseudo-clonal data-set

source lib.sh || exit 1


basedir=data/bamstats
bam=$basedir/bamstats.bam
reffa=$basedir/bamstats.fa
truebamstats=$basedir/bamstats.expected.bamstats

outdir=$(mktemp -d -t $(basename $0).XXXXXX)
outbamstats=$outdir/bamstats.txt
log=$outdir/log.txt

KEEP_TMP=0

cmd="$LOFREQ bamstats -f $reffa  -o $outbamstats $bam"
if ! eval $cmd >> $log 2>&1; then
    echoerror "The following command failed (see $log for more): $cmd"
    exit 1
fi

if ! diff -q $outbamstats $truebamstats; then
    echoerror "Output differs from expected output ($outbamstats differs from $truebamstats)"
    exit 1
else
    echook "Got expected output"
fi



if [ $KEEP_TMP -eq 1 ]; then
    echowarn "Not deleting tmp dir $outdir"
else 
    rm  $outdir/*
    rmdir $outdir
fi

