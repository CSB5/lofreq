#!/bin/bash

# Make sure the parallel wrapper produces the same result as the
# default

source lib.sh || exit 1



BAM=data/dream-icgc-tcga-first10kperchrom-synthetic.challenge.set1.normal.v2.bam
# don't bloody gzip your reference even though samtools happily indexes it
REFFA=data/Homo_sapiens_assembly19.fasta

KEEP_TMP=1
NUM_THREADS=4
DEBUG=0
SIMULATE=0

outdir=$(mktemp -d -t $(basename $0).XXXXXX)
outraw_parallel=$outdir/raw_parallel.vcf
outraw_single=$outdir/raw_single.vcf
log=$outdir/log.txt


cmd="/usr/bin/time -p $LOFREQ call -f $REFFA -o $outraw_single --verbose $BAM"
test $SIMULATE -eq 1 && cmd="echo $cmd"
test $DEBUG -eq 1 && echo "DEBUG: cmd=$cmd" 1>&2
if ! eval $cmd >> $log 2>&1; then
    echoerror "The following command failed (see $log for more): $cmd"
    exit 1
fi

LOFREQ_PARALLEL="$(dirname $LOFREQ)/../lofreq_python/scripts/lofreq2_call_pparallel.py"
cmd="/usr/bin/time -p $LOFREQ_PARALLEL --pp-threads $NUM_THREADS -f $REFFA -o $outraw_parallel --verbose $BAM"
test $SIMULATE -eq 1 && cmd="echo $cmd"
test $DEBUG -eq 1 && echo "DEBUG: cmd=$cmd" 1>&2
if ! eval $cmd >> $log 2>&1; then
    echoerror "The following command failed (see $log for more): $cmd"
    exit 1
fi

if [ $SIMULATE -eq 1 ]; then
    ndiff=0
else
    ndiff=$($LOFREQ vcfset -a complement -1 $outraw_parallel -2 $outraw_single  | grep -c '^[^#]')
fi
if [ $ndiff -ne 0 ]; then
    echoerror "Observed some difference between parallel and single results. Check $outraw_parallel and $outraw_single"
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

