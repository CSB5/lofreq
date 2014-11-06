#!/bin/bash

# FIXME:add-doc

source lib.sh || exit 1

valgrind_log=$(mktemp -t $(basename $0).XXXXXX.valgrind)
vcf_in=data/vcf/CTTGTA_2_remap_razers-i92_peakrem_corr_nodeff.vcf.gz

valgrind --suppression=bgzf_getline.supp --log-file=$valgrind_log --tool=memcheck --leak-check=full $LOFREQ vcfset -a complement -1 $vcf_in -2 $vcf_in >/dev/null || exit 1

test -s $valgrind_log || exit 1

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

