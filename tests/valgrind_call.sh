#!/bin/bash

# FIXME:add-doc

source lib.sh || exit 1

valgrind_log=$(mktemp -t $(basename $0).XXXXXX.valgrind)
vcf_out=$(mktemp -t $(basename $0).XXXXXX.vcf)
rm $vcf_out $valgrind_log

reffa=data/denv2-pseudoclonal/denv2-pseudoclonal_cons.fa
bam=data/denv2-pseudoclonal/denv2-pseudoclonal.bam
# use region, otherwise execution will take ages
region=consensus:100-300 # FIXME
# use hard-coded bonf because:
# 1. can't use auto since we don't have a bed-file
# 2. can't use dynamic since that would eventually call a filtering
# script which valgrind won'tlike
bonf=600

valgrind --log-file=$valgrind_log --tool=memcheck \
    $LOFREQ call -r $region -b $bonf -f $reffa -o $vcf_out $bam || exit 1 

num_snvs=$(grep -cv '^#' $vcf_out)
if [ "$num_snvs" -lt 1 ]; then
    echoerror "Found no SNVs in $vcf_out"
    exit 1
fi

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

#echodebug "FIXME: num_err=$num_err lost_bytes=$lost_bytes num_snvs=$num_snvs"

rm $vcf_out $valgrind_log
