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
bam=../../lofreq-test-data/denv2-simulation/denv2-10haplo.bam
reffa=../../lofreq-test-data/denv2-simulation/denv2-refseq.fa

snv_out_raw=../../lofreq-test-data/denv2-simulation/denv2-10haplo_lofreq-nq-raw.snp
snv_ref=../../lofreq-test-data/denv2-simulation/denv2-10haplo_true-snp.snp

# delete output files from previous run
DEBUG=0
if [ $DEBUG -ne 1 ]; then
    rm -f $snv_out_raw 2>/dev/null
fi
# index bam if necessary
test -s ${bam}.bai || samtools index $bam

# determine bonferroni factors; run LoFreq and filter predictions
bonf=$(lofreq_bonf.py --bam $bam) || exit 1
bonfexp=32169
if [ $bonfexp -ne $bonf ]; then
    echoerror "Expected bonferroni factor to be $bonfexp, but got $bonf. Can't continue" 1>&2
    exit 1
fi
if [ ! -s denv2-10haplo_lofreq-raw.snp ]; then
    lofreq_snpcaller.py --lofreq-q-off --lofreq-nq-on \
        --bonf $bonf -f $reffa -b $bam -o $snv_out_raw || exit 1
else
    echowarn "Reusing snv_out_raw (only useful for debugging)"
fi
echook "Predictions completed."

# test output
nmissing=$(lofreq_diff.py -s $snv_ref -t $snv_out_raw -m uniq_to_1 | wc -l)
nexp=4
if [ $nexp -ne $nmissing ]; then
    echoerror "Number of missing SNVs differs (expected $nexp got $nmissing)"
else
    echook "Got expected number of missing SNVs"
fi
nextra=$(lofreq_diff.py -s $snv_ref -t $snv_out_raw -m uniq_to_2 | wc -l)
nexp=0
if [ $nexp -ne $nextra ]; then
    echoerror "Number of extra SNVs differs (expected $nexp got $nextra)"
else
    echook "Got expected number of extra SNVs"
fi
ncommon=$(lofreq_diff.py -s $snv_ref -t $snv_out_raw -m common | wc -l)
nexp=96
if [ $nexp -ne $ncommon ]; then
    echoerror "Number of common SNVs differs (expected $nexp got $ncommon)"
else
    echook "Got expected number of common SNVs"
fi
