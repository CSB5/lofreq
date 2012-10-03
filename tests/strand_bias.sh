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


snv_raw=../example-data/Ecoli_K12_MG1655_NC_000913_lofreq-bonf.snp


npassfilter_all=$(lofreq_filter.py --strandbias-holmbonf -i $snv_raw | \
    grep -v consensus-var | wc -l)
nexp=11
if [ $nexp -ne $npassfilter_all ]; then
    echoerror "Expected number of (all) SNVs after strand-bias filter differs" 
else
    echook "Got expected number of (all) SNVs after strand-bias filter"
fi

npassfilter_lfv=$(grep -v consensus-var $snv_raw | \
    lofreq_filter.py --strandbias-holmbonf -i - | wc -l)
nexp=2
if [ $nexp -ne $npassfilter_lfv ]; then
    echoerror "Expected number of (low-freq) SNVs after strand-bias filter differs" 
else
    echook "Got expected number of (low-freq) SNVs after strand-bias filter"
fi
