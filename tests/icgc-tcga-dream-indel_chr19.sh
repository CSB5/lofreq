#!/bin/bash

source lib.sh || exit 1

REF=data/icgc-tcga-dream-support/Homo_sapiens_assembly19.fasta
NORMAL=data/icgc-tcga-dream-indel_chr19/chr19.normal_didq_aq.bam
TUMOR=data/icgc-tcga-dream-indel_chr19/chr19.tumor_didq_aq.bam
BED=data/icgc-tcga-dream-indel_chr19/chr19.bed
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

res=$($EVALUATOR -t $TRUTH -v ${outpref}somatic_final.snvs.vcf -m SNV | awk 'END {print $NF}') || exit
echo $res | awk -F, '{prec=$1; rec=$2; if (prec<0.96 || rec<0.96) {status="ERROR"} else {status="OK"} printf "%s: precision=%f recall=%f\n", status, prec, rec > "/dev/stderr"} END {if (status!="OK") {exit 1}}';

res=$($EVALUATOR -t $TRUTH -v ${outpref}somatic_final.indels.vcf -m INDEL | awk 'END {print $NF}') || exit
echo $res | awk -F, '{prec=$1; rec=$2; if (prec<0.96 || rec<0.96) {status="ERROR"} else {status="OK"} printf "%s: precision=%f recall=%f\n", status, prec, rec > "/dev/stderr"} END {if (status!="OK) {exit 1}"}'



            