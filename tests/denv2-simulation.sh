#!/bin/bash

# Call SNVs on simulated data and make sure we got the expected number
# of SNVs

source lib.sh || exit 1

basedir=data/denv2-simulation
bam=$basedir/denv2-10haplo.bam
reffa=$basedir/denv2-refseq.fa
truesnv=$basedir/denv2-10haplo_true-snp.vcf.gz
# samtools mpileup $bam | wc -l;# *3
bonf=32169

outdir=$(mktemp -d -t $(basename $0).XXXXXX)
outraw_def=$outdir/raw_def.vcf
outfinal_def=$outdir/final_def.vcf
outraw_nomq=$outdir/raw_nomq.vcf
outfinal_nomq=$outdir/final_nomq.vcf
log=$outdir/log.txt

KEEP_TMP=0
if [ $KEEP_TMP -eq 1 ]; then
    echowarn "Keeping tmp dir $outdir"
fi

cmd="$LOFREQ call -B -b $bonf -f $reffa -o $outraw_def $bam"
if ! eval $cmd >> $log 2>&1; then
    echoerror "The following command failed (see $log for more): $cmd"
    exit 1
fi
cmd="$LOFREQ filter --only-passed -i $outraw_def -o $outfinal_def"
if ! eval $cmd >> $log 2>&1; then
    echoerror "The following command failed (see $log for more): $cmd"
    exit 1
fi


cmd="$LOFREQ call -B -b $bonf -f $reffa -o $outraw_nomq -J $bam"
if ! eval $cmd >> $log 2>&1; then
    echoerror "The following command failed (see $log for more): $cmd"
    exit 1
fi
cmd="$LOFREQ filter --only-passed -i $outraw_nomq -o $outfinal_nomq"
if ! eval $cmd >> $log 2>&1; then
    echoerror "The following command failed (see $log for more): $cmd"
    exit 1
fi


#nexp=$(grep -v -c '^#' $truesnv)
#nfinal_def=$(grep -v -c '^#' $outfinal_def)
#nfinal_nomq=$(grep -v -c '^#' $outfinal_nomq)
#echodebug "nexp=$nexp nfinal_def=$nfinal_def $nfinal_nomq=$nfinal_nomq"


ndiff=$($LOFREQ vcfset -a complement --only-passed -1 $outfinal_def -2 $truesnv  | grep -c '^[^#]')
if [ $ndiff -ne 0 ]; then
    echoerror "Found extra SNVs in default predictions, which are not part of the list of true SNVs"
    exit 1
fi
ndiff=$($LOFREQ vcfset -a complement --only-passed -2 $outfinal_def -1 $truesnv  | grep -c '^[^#]')
nexp=15
# BAQ on: 19
# BAQ off: 15
if [ $ndiff -ne $nexp ]; then
    echoerror "Expected $nexp missing SNVs in default predictions but got $ndiff"
    exit 1
fi



ndiff=$($LOFREQ vcfset -a complement --only-passed -1 $outfinal_nomq -2 $truesnv  | grep -c '^[^#]')
if [ $ndiff -ne 0 ]; then
    echoerror "Found extra SNVs in no-mq predictions, which are not part of the list of true SNVs"
    exit 1
fi
ndiff=$($LOFREQ vcfset -a complement --only-passed -2 $outfinal_nomq -1 $truesnv  | grep -c '^[^#]')
nexp=11
# BAQ on: 14
# BAQ off: 11
if [ $ndiff -ne $nexp ]; then
    echoerror "Expected $nexp missing SNVs in no-mq predictions but got $ndiff"
    exit 1
fi


# FIXME outfinal should not look different, i.e. filtering shouldn't do much/anything.
# see /home/wilma/snpcaller/lofreq/lofreq-sourceforge.git/tests/denv2-simulation.sh 

echook "Tests passed"

if [ $KEEP_TMP -ne 1 ]; then
    rm  $outdir/*
    rmdir $outdir
fi

