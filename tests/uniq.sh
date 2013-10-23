#!/bin/bash

# FIXME:add-doc

source lib.sh || exit 1

# test vs self prediction should give zero results

bam=data/denv2-simulation/denv2-10haplo.bam
vcf_in=data/denv2-simulation/denv2-10haplo_true-snp.vcf.gz
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


vcf_in=data/vcf/denv2-10haplo-fake-filter-only-and-indels.vcf.gz
$LOFREQ uniq -v $vcf_in $bam -o $vcf_out || exit 1
num_snvs=$(grep -cv '^#' $vcf_out)
if [ "$num_snvs" -ne 0 ]; then
    echoerror "Expected zero SNVs when checking against indels and filtered variants only but got $num_snvs"
    exit 1
else
    echook "Got zero SNVs when checking indels and filtered variants only"
fi
rm $vcf_out


vcf_in=data/somatic/hg19_chr22_true_snv.vcf.gz
bam=data/somatic/CHH966-tumor-100x-10pur-hg19.chr22-bed-only.bam
$LOFREQ uniq -v $vcf_in $bam -o $vcf_out || exit 1
# previously 4, but now 2 true snvs in vcf_in, which both should be unique
num_snvs=$(grep -cv '^#' $vcf_out)
if [ "$num_snvs" -ne 2 ]; then
    echoerror "Expected two SNVs from somatic check but got $num_snvs"
    exit 1
else
    echook "Got expected number of SNVs from somatic check"
fi

rm $vcf_out
#echo $vcf_out


