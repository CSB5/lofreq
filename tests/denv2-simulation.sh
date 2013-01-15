#!/bin/bash

source lib.sh || exit 1

bam=../../lofreq-test-data/denv2-simulation/denv2-10haplo.bam
reffa=../../lofreq-test-data/denv2-simulation/denv2-refseq.fa

snv_out_raw=../../lofreq-test-data/denv2-simulation/denv2-10haplo_lofreq-raw.snp
snv_out_nonjoined_raw=../../lofreq-test-data/denv2-simulation/denv2-10haplo_lofreq-nonjoined-raw.snp
snv_out_sbf=../../lofreq-test-data/denv2-simulation/denv2-10haplo_lofreq-sbf.snp
snv_ref=../../lofreq-test-data/denv2-simulation/denv2-10haplo_true-snp.snp

# delete output files from previous run
DEBUG=0
if [ $DEBUG -ne 1 ]; then
    rm -f $snv_out_raw $snv_out_nonjoined_raw $snv_out_sbf 2>/dev/null
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
if [ ! -s $snv_out_raw ]; then
    lofreq_snpcaller.py --bonf $bonf -f $reffa -b $bam -o $snv_out_raw || exit 1
else
    echowarn "Reusing snv_out_raw (only useful for debugging)"
fi

if [ ! -s $snv_out_nonjoined_raw ]; then
    lofreq_snpcaller.py --dont-join-mapq-and-baseq --bonf $bonf -f $reffa -b $bam -o $snv_out_nonjoined_raw || exit 1
else
    echowarn "Reusing snv_out_nonjoined_raw (only useful for debugging)"
fi

nnonjoined=$(cat $snv_out_nonjoined_raw | wc -l)
norig=$(cat $snv_out_raw | wc -l)
if [ $nnonjoined -lt $norig ]; then
    echoerror "Was expecting more or equal number of SNVs when using nonjoined-baseq-mapq"
elif [ $nnonjoined -eq 0 ]; then
    echoerror "Got 0 SNVs when using nonjoined-baseq-mapq"
else
    echook "Got Nonjoined-BaseQ-MapQ results look ok"
fi


if [ ! -s $snv_out_sbf ]; then
    lofreq_filter.py --strandbias-holmbonf \
        -i $snv_out_raw \
        -o $snv_out_sbf || exit 1
else
    echowarn "Reusing $snv_out_sbf (only useful for debugging)" 1>&2
fi
echook "Predictions completed."

# test non-joined output
nmissing=$(lofreq_diff.py -s $snv_ref -t $snv_out_nonjoined_raw -m uniq_to_1 | wc -l)
nexp=14
if [ $nexp -ne $nmissing ]; then
    echoerror "Number of missing SNVs differs (expected $nexp got $nmissing)"
else
    echook "Got expected number of missing SNVs"
fi
nextra=$(lofreq_diff.py -s $snv_ref -t $snv_out_nonjoined_raw -m uniq_to_2 | wc -l)
nexp=0
if [ $nexp -ne $nextra ]; then
    echoerror "Number of extra SNVs differs (expected $nexp got $nextra)"
else
    echook "Got expected number of extra SNVs"
fi
ncommon=$(lofreq_diff.py -s $snv_ref -t $snv_out_nonjoined_raw -m common | wc -l)
nexp=86
if [ $nexp -ne $ncommon ]; then
    echoerror "Number of common SNVs differs (expected $nexp got $ncommon)"
else
    echook "Got expected number of common SNVs"
fi


# test normal output
nmissing=$(lofreq_diff.py -s $snv_ref -t $snv_out_raw -m uniq_to_1 | wc -l)
nexp=19
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
nexp=81
if [ $nexp -ne $ncommon ]; then
    echoerror "Number of common SNVs differs (expected $nexp got $ncommon)"
else
    echook "Got expected number of common SNVs"
fi


# SBF filtering doesn't do anything here since the simulation doesn't
# introduce strandbias lofreq_filter does correction and changes info.
# Only look at first four fields therefore
md5raw=$(cut -f1-4 -d ' ' $snv_out_raw | $md5)
md5sbf=$(cut -f1-4 -d ' ' $snv_out_sbf | $md5)
if [ "$md5raw" != "$md5raw" ]; then
    echoerror "SB-filtering changed number of SNVs"
else
    echook "SB-filtering didn't do much as expected here (simulation)"
fi
