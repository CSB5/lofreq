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
bam=../example-data/denv2-pseudoclonal/denv2-pseudoclonal.bam
reffa=../example-data/denv2-pseudoclonal/denv2-pseudoclonal_cons.fa
bed=../example-data/denv2-pseudoclonal/denv2-pseudoclonal_incl.bed

snv_out_raw=../example-data/denv2-pseudoclonal/denv2-pseudoclonal_lofreq-raw.snp
snv_out_sbf=../example-data/denv2-pseudoclonal/denv2-pseudoclonal_lofreq-sbf.snp
snv_ref=../example-data/denv2-pseudoclonal/denv2-pseudoclonal_true-snp.snp

# delete output files from previous run
DEBUG=0
if [ $DEBUG -ne 1 ]; then
    rm -f $snv_out_raw $snv_out_sbf 2>/dev/null
fi

# index bam if necessary
test -s ${bam}.bai || samtools index $bam

# determine bonferroni factors; run LoFreq and filter predictions
bonf=$(lofreq_bonf.py -i $bam -b $bed) || exit 1
bonfexp=29727
if [ $bonfexp -ne $bonf ]; then
    echoerror "Expected bonferroni factor to be $bonfexp, but got $bonf. Can't continue" 1>&2
    exit 1
fi
if [ ! -s $snv_out_raw ]; then
    samtools mpileup -d 100000 -E -f $reffa -l $bed $bam | \
        lofreq_snpcaller.py -i - -b $bonf \
        -o $snv_out_raw || exit 1
else
    echowarn "Reusing $snv_out_raw (only useful for debugging)"
fi
if [ ! -s $snv_out_sbf ]; then
    lofreq_filter.py --strandbias-holmbonf \
        -i $snv_out_raw \
        -o $snv_out_sbf || exit 1
else
    echowarn "Reusing $snv_out_sbf (only useful for debugging)"
fi
echook "Predictions completed."

# test raw output
nmissing=$(lofreq_diff.py -s $snv_ref -t $snv_out_raw -m uniq_to_1 | wc -l)
nexp=58
if [ $nexp -ne $nmissing ]; then
    echoerror "Number of missing SNVs differs (expected $nexp got $nmissing)"
else
    echook "Got expected number of missing SNVs"
fi
nextra=$(lofreq_diff.py -s $snv_ref -t $snv_out_raw -m uniq_to_2 | wc -l)
#nexp=18
nexp=0
if [ $nexp -ne $nextra ]; then
    echoerror "Number of extra SNVs differs (expected $nexp got $nextra)"
else
    echook "Got expected number of extra SNVs"
fi
ncommon=$(lofreq_diff.py -s $snv_ref -t $snv_out_raw -m common | wc -l)
nexp=229
if [ $nexp -ne $ncommon ]; then
    echoerror "Number of common SNVs differs (expected $nexp got $ncommon)"
else
    echook "Got expected number of common SNVs"
fi


# SBF filtering doesn't do anything here since the simulation doesn't
# introduce strandbias lofreq_filter does correction and changes info.
# Only look at first four fields therefore
md5raw=$(cut -f1-4 -d ' ' $snv_out_raw | md5sum)
md5sbf=$(cut -f1-4 -d ' ' $snv_out_sbf | md5sum)
if [ "$md5raw" != "$md5raw" ]; then
    echoerror "SB-filtering changed number of SNVs"
else
    echook "SB-filtering didn't do much as expected here (explicitely removed)"
fi
