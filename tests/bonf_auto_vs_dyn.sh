#!/bin/bash

# Call SNVs on a BAM fule with full coverage and change bonf settings.
# Different settings (auto, dynamic and hard-coded) should give
# identical results here.

source lib.sh || exit 1

basedir=data/denv2-pseudoclonal
bed=$basedir/denv2-pseudoclonal_incl.bed
bam=$basedir/denv2-pseudoclonal.bam
reffa=$basedir/denv2-pseudoclonal_cons.fa

outdir=$(mktemp -d -t $(basename $0).XXXXXX)
out_auto=$outdir/snv_auto.vcf
out_dynamic=$outdir/snv_dynamic.vcf
# bed_len.sh $be;# = 9909 * 3 = 29727
out_29727=$outdir/snv_29727.vcf
log=$outdir/log.txt

KEEP_TMP=0

cmd="$LOFREQ call -l $bed -b auto -f $reffa -o $out_auto $bam"
if ! eval $cmd >> $log 2>&1; then
    echoerror "The following command failed (see $log for more): $cmd"
    exit 1
fi

cmd="$LOFREQ call -l $bed -b dynamic -f $reffa -o $out_dynamic $bam"
if ! eval $cmd >> $log 2>&1; then
    echoerror "The following command failed (see $log for more): $cmd"
    exit 1
fi

cmd="$LOFREQ call -l $bed -b 29727 -f $reffa -o $out_29727 $bam"
if ! eval $cmd >> $log 2>&1; then
    echoerror "The following command failed (see $log for more): $cmd"
    exit 1
fi

#echodebug "All calls done. No checking results"

# make sure we got at least some snvs
# 
if [ $(grep -c '^[^#]' $out_auto) -eq 0 ]; then
    echoerror "No SNVs predicted"
    exit 1
fi

ndiff=$(lofreq2_vcfset.py -a complement --ign-filtered -1 $out_auto -2 $out_dynamic  | grep -c '^[^#]') || exit 1
if [ $ndiff -ne 0 ]; then
    echoerror "Found differences between bonf auto and bonf dynamic outputs"
    exit 1
fi
ndiff=$(lofreq2_vcfset.py -a complement --ign-filtered -2 $out_dynamic -1 $out_auto  | grep -c '^[^#]') || exit 1
if [ $ndiff -ne 0 ]; then
    echoerror "Found differences between bonf auto and bonf dynamic outputs"
    exit 1
fi

ndiff=$(lofreq2_vcfset.py -a complement --ign-filtered -1 $out_auto -2 $out_29727  | grep -c '^[^#]') || exit 1
if [ $ndiff -ne 0 ]; then
    echoerror "Found differences between bonf auto and bonf 29727 outputs"
    exit 1
fi

ndiff=$(lofreq2_vcfset.py -a complement --ign-filtered -2 $out_29727 -1 $out_auto  | grep -c '^[^#]') || exit 1
if [ $ndiff -ne 0 ]; then
    echoerror "Found differences between bonf auto and bonf 29727 outputs"
    exit 1
fi

echook "Tests passed"

if [ $KEEP_TMP -eq 1 ]; then
    echowarn "Not deleting tmp dir $outdir"
else 
    rm  $outdir/*
    rmdir $outdir
fi
