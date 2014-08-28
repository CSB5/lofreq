#!/bin/bash

source lib.sh || exit 1

KEEP_TMP=0
REF=data/icgc-tcga-dream-support/Homo_sapiens_assembly19.fasta
NORMAL=data/icgc-tcga-dream-indel_chr19/chr19.normal_didq_aq.bam
TUMOR=data/icgc-tcga-dream-indel_chr19/chr19.tumor_didq_aq.bam
BED=data/icgc-tcga-dream-indel_chr19/chr19.bed
#BED=data/icgc-tcga-dream-indel_chr19/chr19-debug.bed
EVALUATOR=data/icgc-tcga-dream-support/evaluator.py
TRUTH=data/icgc-tcga-dream-indel_chr19/chr19.truth.vcf.gz

outdir=$(mktemp -d -t $(basename $0).XXXXXX)
outpref=$outdir/lofreq_test
log=$outdir/log.txt

cmd="$LOFREQ somatic -f $REF --threads $threads -n $NORMAL -t $TUMOR -o $outpref -l $BED --verbose"
echodebug "cmd=$cmd"
if ! eval $cmd > $log 2>&1; then
    echoerror "LoFreq failed. Check logfile $log. Command was $cmd"
    exit 1
fi

num_err=0

res_ll=$($EVALUATOR -t $TRUTH -v ${outpref}somatic_final.snvs.vcf -m SNV | awk 'END {print $NF}') || exit 1
res=$(echo $res_ll | \
  awk -F, '{prec=$1; rec=$2; if (prec<0.96 || rec<0.96) {status="ERROR"} else {status="OK"} printf "%s: precision=%f recall=%f\n", status, prec, rec}') || exit 1
if echo $res | grep -q ERROR; then
   let num_err=num_err+1
fi
echo $res 1>&2

res_ll=$($EVALUATOR -t $TRUTH -v ${outpref}somatic_final.indels.vcf -m INDEL | awk 'END {print $NF}') || exit 1
res=$(echo $res_ll | \
  awk -F, '{prec=$1; rec=$2; if (prec<0.90 || rec<0.50) {status="ERROR"} else {status="OK"} printf "%s: precision=%f recall=%f\n", status, prec, rec}') || exit 1
if echo $res | grep -q ERROR; then
   let num_err=num_err+1
fi
echo $res 1>&2

if [ $KEEP_TMP -ne 1 ] && [ $num_err -eq 0 ]; then
    test -d $outdir && rm -rf $outdir
else
    echowarn "Not deleting temporary output directory $outdir"
fi
if [ $num_err -ne 0 ]; then
    exit 1
fi
