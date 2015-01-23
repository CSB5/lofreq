#!/bin/bash
source lib.sh || exit 1

KEEP_TMP=0
BASE=data/icgc-tcga-dream-testproject/
#BASE=/projects/wilma/SOMATIC/dream-challenge/testproject/
#BASE=/mnt/userArchive/wilma/projects/somatic/testproject/
REF=data/icgc-tcga-dream-support/Homo_sapiens_assembly19.fasta
TUMOR=${BASE}/tumor.chr20.bam
NORMAL=${BASE}/normal.chr20.bam
TRUTH=${BASE}/truth.chr20.vcf.gz
BED=${BASE}/chr20.bed
#EVALUATOR=/projects/wilma/SOMATIC/dream-challenge/tools/bamsurgeon.git/etc/evaluator.py
EVALUATOR=data/icgc-tcga-dream-support/evaluator.py
DEBUG=0

# threads=16; echoinfo "overwriting default threads to $threads"


for f in $REF $TUMOR $NORMAL $TRUTH $EVALUATOR; do
    if [ ! -s $f ]; then
        echoerror "Essential file $f missing"
        exit 1
    fi
done
out_pref=$(mktemp -t $(basename $0).XXXXXX)
log=${out_pref}.exec.log
vcf_out=${out_pref}somatic_final.snvs.vcf.gz
if [ $DEBUG -eq 1 ]; then
    cp ${BASE}/snvs/lofreq/beta-4-8-g7b8b334-dirty_somatic_final.vcf $vcf_out
else
    cmd="$LOFREQ somatic -l $BED -n $NORMAL -t $TUMOR -f $REF -o $out_pref --threads $threads"
    if ! eval $cmd > $log 2>&1; then
        echoerror "LoFreq failed. Check log $log and files with prefix $out_pref"
        exit 1
    fi
    echoinfo "lofreq somatic run completed. now checking results"
fi

num_err=0
# use bamsurgeon evaluator
#
# example output
# alterantive to using evaluator is to run lofreq vcfset on a truth file only containing SNVs
# tpcount, fpcount, subrecs, trurecs:
# 1389 15 1404 1445
# precision, recall, F1 score: 0.989316239316,0.96124567474,0.975078975079
title="snvs before dbsnp removal"
res_ll=$($EVALUATOR -t $TRUTH -v $vcf_out -m SNV | awk 'END {print $NF}') || exit

res=$(echo $res_ll | \
    awk -F, '{prec=$1; rec=$2; if (prec<0.98 || rec<0.96) {status="ERROR"} else {status="OK"} printf "%s: precision=%f recall=%f\n", status, prec, rec}') || exit 1
if echo $res | grep -q ERROR; then
   let num_err=num_err+1
fi
echo "$title: "$res 1>&2


if [ $KEEP_TMP -ne 1 ] && [ $num_err -eq 0 ]; then
    test -d $outdir && rm -rf $outdir
else
    echowarn "Not deleting temporary output directory $outdir"
fi
if [ $num_err -ne 0 ]; then
    exit 1
fi

