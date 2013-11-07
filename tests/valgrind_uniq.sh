#!/bin/bash

# FIXME:add-doc

source lib.sh || exit 1

valgrind_log=$(mktemp -t $(basename $0).XXXXXX.valgrind)
vcf_out=$(mktemp -t $(basename $0).XXXXXX.vcf)
rm $vcf_out $valgrind_log

# FIXME better to use somatic SNVs
bam=data/denv2-simulation/denv2-10haplo.bam
vcf=data/denv2-simulation/denv2-10haplo_true-snp.vcf.gz

# use only head. otherwise too slow
zcat $vcf | head | valgrind  --log-file=$valgrind_log --tool=memcheck \
    $LOFREQ uniq -v - $bam -o $vcf_out || exit 1
 

num_err=$(grep 'ERROR SUMMARY' $valgrind_log | grep -cv ': 0 errors')
if [ "$num_err" -ne 0 ]; then
    echoerror "Found errors in Valgrind output $valgrind_log"
    exit 1
else
    echook "No errors found in Valgrind output"
fi

lost_bytes=$(grep 'lost' $valgrind_log | grep -cv ': 0 bytes in 0 blocks')
if [ "$lost_bytes" -ne 0 ]; then
    echoerror "Found lost bytes in Valgrind output $valgrind_log" || exit 1
    exit 1
else
    echook "No lost bytes found in Valgrind output"
fi



rm $vcf_out $valgrind_log
