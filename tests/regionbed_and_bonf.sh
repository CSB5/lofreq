#!/bin/bash

echoerror() {
    echo "ERROR: $1" 1>&2
}
echook() {
    echo "OK: $1" 1>&2
}
echowarn() {
    echo "WARN: $1" 1>&2
}


bam=../../lofreq-test-data/gastric-cancer-wgs/WHG001-final.bam
bed=${bam%.bam}.bed

rm $bed 2>/dev/null

lofreq_regionbed.py -i $bam -o $bed || exit 1
bonf=$(lofreq_bonf.py --bed $bed) || exit 1 

exp=9287610339
if [ $exp -ne $bonf ]; then
    echoerror "Expected Bonferroni factor differs" 
else
    echook "Got expected Bonferroni factor"
fi

