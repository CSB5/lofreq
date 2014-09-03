#!/bin/bash
source lib.sh || exit 1

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

for f in $REF $TUMOR $NORMAL $TRUTH $EVALUATOR; do
    if [ ! -s $f ]; then
        echoerror "Essential file $f missing"
        exit 1
    fi
done
out_pref=$(mktemp -t $(basename $0).XXXXXX)
vcf_out=${out_pref}somatic_final.snvs.vcf
if [ $DEBUG -eq 1 ]; then
    cp ${BASE}/snvs/lofreq/beta-4-8-g7b8b334-dirty_somatic_final.vcf $vcf_out
else
    $LOFREQ somatic -l $BED -n $NORMAL -t $TUMOR -f $REF -o $out_pref --threads $threads || exit 1
    echoinfo "lofreq somatic run completed. now checking results"
fi

# use bamsurgeon evaluator
#
# example output
# alterantive to using evaluator is to run lofreq vcfset on a truth file only containing SNVs
# tpcount, fpcount, subrecs, trurecs:
# 1389 15 1404 1445
# precision, recall, F1 score: 0.989316239316,0.96124567474,0.975078975079
res=$($EVALUATOR -t $TRUTH -v $vcf_out -m SNV | awk 'END {print $NF}') || exit
echo $res | awk -F, '{prec=$1; rec=$2; if (prec<0.96 || rec<0.96) {status="ERROR"} else {status="OK"} printf "%s: precision=%f recall=%f\n", status, prec, rec > "/dev/stderr"}'

#rm $vcf_out

