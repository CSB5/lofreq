#!/bin/bash
source lib.sh || exit 1

REF=/mnt/userArchive/wilma/projects/somatic/testproject/refseq/Homo_sapiens_assembly19.fasta
TUMOR=/mnt/userArchive/wilma/projects/somatic/testproject/tumor.chr20.bam
NORMAL=/mnt/userArchive/wilma/projects/somatic/testproject/normal.chr20.bam
TRUTH=/mnt/userArchive/wilma/projects/somatic/testproject/truth.chr20.vcf.gz
BED=/mnt/userArchive/wilma/projects/somatic/testproject/chr20.bed
EVALUATOR=/projects/wilma/SOMATIC/dream-challenge/tools/bamsurgeon.git/etc/evaluator.py
DEBUG=0

for f in $REF $TUMOR $NORMAL $TRUTH $EVLUATOR; do
    if [ ! -s $f ]; then
        echoerror "Essential file $f missing"
        exit 1
    fi
done
vcf_out=$(mktemp -t $(basename $0).XXXXXX.vcf)
if [ $DEBUG -eq 1 ]; then
    cp /mnt/userArchive/wilma/projects/somatic/testproject/snvs/lofreq/beta-4-8-g7b8b334-dirty_somatic_final.vcf $vcf_out
else
    $LOFREQ somatic -l $BED -n $NORMAL -t $TUMOR -f $REF -o $vcf_out --threads $threads || exit 1
    echoinfo "lofreq somatic run completed. now checking results"
fi

# use bamsurgeon evaluator
#
# example output
# alterantive to using evaluator is to run lofreq vcfset on a truth file only containing SNVs
# tpcount, fpcount, subrecs, trurecs:
# 1389 15 1404 1445
# precision, recall, F1 score: 0.989316239316,0.96124567474,0.975078975079
res=$($EVALUATOR -t $TRUTH -v $vcf_out -m SNV | awk 'END {print $NF}')
echo $res | awk -F, '{prec=$1; rec=$2; if (prec<0.96 || rec<0.96) {status="ERROR"} else {status="OK"} printf "%s: precision=%f recall=%f\n", status, prec, rec > "/dev/stderr"}'

#rm $vcf_out

