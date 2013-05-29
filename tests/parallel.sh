#!/bin/bash

# Make sure the parallel wrapper produces the same result as the
# default

source lib.sh || exit 1


basedir=data/denv2-pseudoclonal
bam=$basedir/denv2-pseudoclonal.bam
reffa=$basedir/denv2-pseudoclonal_cons.fa
bed=$basedir/denv2-pseudoclonal_incl.bed
#truesnv=$basedir/denv2-pseudoclonal_true-snp.vcf

basedir=data/somatic
bam=$basedir/CHH966-tumor-100x-10pur-hg19.chr22-bed-only.bam
reffa=$basedir/hg19_chr22.fa
bed=$basedir/SeqCap_EZ_Exome_v3_primary_lib_extend_no_overlap_minus300.chr22.bed


outdir=$(mktemp -d -t $(basename $0).XXXXXX)
outraw_parallel=$outdir/raw_parallel.vcf
outraw_single=$outdir/raw_single.vcf
log=$outdir/log.txt

KEEP_TMP=0
NUM_THREADS=4

cmd="$LOFREQ call -f $reffa -l $bed -b dynamic -o $outraw_single --verbose $bam"
if ! eval $cmd >> $log 2>&1; then
    echoerror "The following command failed (see $log for more): $cmd"
    exit 1
fi

cmd="$(dirname $LOFREQ)/../lofreq_python/scripts/lofreq2_call_parallel.py -n $NUM_THREADS -f $reffa -l $bed -b dynamic -o $outraw_parallel --verbose $bam"
if ! eval $cmd >> $log 2>&1; then
    echoerror "The following command failed (see $log for more): $cmd"
    exit 1
fi


ndiff=$($LOFREQ vcfset -a complement -1 $outraw_parallel -2 $outraw_single  | grep -c '^[^#]')
if [ $ndiff -ne 0 ]; then
    echoerror "Observed some difference between parallel and single results. Check $outraw_parallel and $outraw_single"
    exit 1
else
    echook "Parallel and single run give identical results."
    echodebug "EXIT without deleting $outdir"; exit 1
fi



if [ $KEEP_TMP -eq 1 ]; then
    echowarn "Not deleting tmp dir $outdir"
else 
    rm  $outdir/*
    rmdir $outdir
fi

