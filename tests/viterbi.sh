#!/bin/bash

source lib.sh || exit 1

BASEDIR=data/viterbi/
REF=$BASEDIR/NC_011770.fa
BAM=$BASEDIR/pseudomonas_pair_screwed_up_cigar.bam


# input contains two reads with near random cigar strings
# that are in fact perfect matches

ncorr=$($LOFREQ viterbi -f $REF $BAM | samtools view - 2>/dev/null | grep -cw 75M) || exit 1
if [ $ncorr != "2" ]; then
    echoerror "Expected two fixed input reads but got $ncorr"
    exit 1
else
    echook "All reads correctly realigned"
fi
