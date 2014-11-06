#!/bin/bash

# FIXME:add-doc

source lib.sh || exit 1

valgrind_log=$(mktemp -t $(basename $0).XXXXXX.valgrind)
vcf_out=$(mktemp -t $(basename $0).XXXXXX.vcf)
rm $vcf_out $valgrind_log

reffa=data/denv2-pseudoclonal/denv2-pseudoclonal_cons.fa
bam=data/denv2-pseudoclonal/denv2-pseudoclonal.bam
#reffa=./data/somatic/hg19_chr22.fa.gz
#bam=./data/somatic/CHH966-tumor-100x-10pur-hg19.chr22-bed-only.bam
# use region, otherwise execution will take ages
#region_arg=consensus:100-300 # FIXME

# use hard-coded bonf because:
# 1. can't use auto since we don't have a bed-file
# 2. can't use dynamic since that would eventually call a filtering
# script which valgrind won'tlike
bonf_arg="-b 308334"
#region_arg="-r consensus:1-1000"

valgrind --suppressions=faidx_fetch_seq.supp --log-file=$valgrind_log --tool=memcheck --leak-check=full \
    $LOFREQ call $region_arg $bonf_arg --no-default-filter -f $reffa -o $vcf_out $bam || exit 1

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



rm $vcf_out $valgrind_log
