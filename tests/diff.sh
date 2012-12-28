#!/bin/bash

source lib.sh || exit 1

basedir=../../lofreq-test-data//alnoffset/
snp1=$basedir/CTTGTA_2_remap_razers-i92_peakrem_corr.snp
snp2=$basedir/CTTGTA_2_remap_razers-i92_peakrem_corr_sbhbf.snp

test -s $snp1 || exit 1

# to self should be zero
ndiffs=$(lofreq_diff.py -s $snp1 -t $snp1 -m uniq_to_2 | wc -l)
nexp=0
if [ $nexp -ne $ndiffs ]; then
    echoerror "Number of SNV differences should be zero when compared against self"
else
    echook "Got expected number of diff SNVs"
fi

# some should be "uniq" in 1 if compared against sb-filtered
ndiffs=$(lofreq_diff.py -s $snp1 -t $snp2 -m uniq_to_1 | wc -l)
nexp=11
if [ $nexp -ne $ndiffs ]; then
    echoerror "Number of expected SNV differences differs"
else
    echook "Got expected number of diff SNVs"
fi

# none should be uniq in sb-filtered
ndiffs=$(lofreq_diff.py -s $snp1 -t $snp2 -m uniq_to_2 | wc -l)
nexp=0      
if [ $nexp -ne $ndiffs ]; then
    echoerror "Filtered SNV file should not have 'uniq' SNVs"
else
    echook "Got expected number of diff SNVs"
fi
