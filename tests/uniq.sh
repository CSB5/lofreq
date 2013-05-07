#!/bin/bash

# FIXME:add-doc

source lib.sh || exit 1


# tst vs self prediction should give zero results

bam=data/denv2-simulation/denv2-10haplo.bam
vcf_in=data/denv2-simulation/denv2-10haplo_true-snp.vcf
vcf_out=$(mktemp -t $(basename $0).XXXXXX.vcf)
rm $vcf_out

$LOFREQ uniq -v $vcf_in $bam -o $vcf_out || exit 1
num_snvs=$(grep -cv '^#' $vcf_out)
if [ "$num_snvs" -ne 0 ]; then
    echoerror "Expected zero SNVs when checking variants predicted from same BAM but got $num_snvs"
    exit 1
else
    echook "Got zero SNVs during self-comparison, as expected"
fi
rm $vcf_out


vcf_in=data/vcf/denv2-10haplo-fake-filter-only-and-indels.vcf
$LOFREQ uniq -v $vcf_in $bam -o $vcf_out || exit 1
num_snvs=$(grep -cv '^#' $vcf_out)
if [ "$num_snvs" -ne 0 ]; then
    echoerror "Expected zero SNVs when checking against indels and filtered variants only but got $num_snvs"
    exit 1
else
    echook "Got zero SNVs when checking indels and filtered variants only"
fi
rm $vcf_out
