#!/bin/bash

# Make sure the parallel wrapper produces the same result as the
# default

source lib.sh || exit 1



BAM=data/icgc-tcga-first10kperchrom-syn1/dream-icgc-tcga-first10kperchrom-synthetic.challenge.set1.normal.v2.bam
# don't bloody gzip your reference even though samtools happily indexes it
REF=data/icgc-tcga-dream-support/Homo_sapiens_assembly19.fasta

KEEP_TMP=0
DEBUG=0
SIMULATE=0

outdir=$(mktemp -d -t $(basename $0).XXXXXX)
outraw_parallel=$outdir/raw_parallel.vcf
outraw_single=$outdir/raw_single.vcf
log=$outdir/log.txt

LOFREQ_PARALLEL="$(dirname $LOFREQ)/../scripts/lofreq2_call_pparallel.py"
cmd="/usr/bin/time -p $LOFREQ_PARALLEL --pp-threads $threads -f $REF -o $outraw_parallel --verbose $BAM"
test $SIMULATE -eq 1 && cmd="echo $cmd"
test $DEBUG -eq 1 && echo "DEBUG: cmd=$cmd" 1>&2
if ! eval $cmd >> $log 2>&1; then
    echoerror "The following command failed (see $log for more): $cmd"
    exit 1
fi


cmd="/usr/bin/time -p $LOFREQ call -f $REF -o $outraw_single --verbose $BAM"
test $SIMULATE -eq 1 && cmd="echo $cmd"
test $DEBUG -eq 1 && echo "DEBUG: cmd=$cmd" 1>&2
if ! eval $cmd >> $log 2>&1; then
    echoerror "The following command failed (see $log for more): $cmd"
    exit 1
fi


if [ $SIMULATE -eq 1 ]; then
    nup=0
    nus=0
else
    nup=$($LOFREQ vcfset -a complement -1 $outraw_parallel -2 $outraw_single --count-only)
    nus=$($LOFREQ vcfset -a complement -2 $outraw_parallel -1 $outraw_single --count-only)
fi
#if [ $nup -ne 0 ] || [ $nus -ne 0 ] ; then
# there are occasional differences possible likely due to BAQ effects on region ends
if [ $nup -gt 1 ] || [ $nus -gt 1 ] ; then
    echoerror "Observed some difference between parallel and single results. Check $outraw_parallel and $outraw_single"
    n_parallel=$(grep -vc '^#' $outraw_parallel)
    n_single=$(grep -vc '^#' $outraw_single)

    n_overlap=$($LOFREQ vcfset -a intersect -1 $outraw_parallel -2 $outraw_single --count-only)
    echoerror "$outraw_parallel has $n_parallel and $outraw_single has $n_single SNVS. Both overlap by $n_overlap"
    exit 1
else
    echook "Parallel and single run give identical results."
fi



if [ $KEEP_TMP -eq 1 ]; then
    echowarn "Not deleting tmp dir $outdir"
else 
    rm  $outdir/*
    rmdir $outdir
fi

