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




# no indels!
vcf_in=data/vcf/CTTGTA_2_remap_razers-i92_peakrem_corr_nodeff.vcf.gz
bam=data/denv2-dpcr-validated/GGCTAC_2_remap_razers-i92_peakrem_corr.bam

# in == out with detlim
#
num_in=$(zgrep -cv '^#' $vcf_in)
cmd="$LOFREQ uniq -v $vcf_in $bam --use-det-lim -o -"
num_out=$(eval $cmd | grep -vc '^#') || exit 1
if [ "$num_in" -ne "$num_out" ]; then
    echoerror "Expected same number of in and output vars when using --use-det-lim but go $num_in and $num_out resp. (cmd was $cmd)"
fi

# UQ= present even with --output-all
cmd="$LOFREQ uniq -v $vcf_in $bam --output-all -o -"
eval $cmd | grep -q 'UQ=' || echoerror "No UQ markup found"

# in gt out in default mode
num_in=$(zgrep -cv '^#' $vcf_in)
cmd="$LOFREQ uniq -v $vcf_in $bam -o -"
num_out=$(eval $cmd | grep -vc '^#') || exit 1
if [ "$num_in" -le "$num_out" ]; then
    echoerror "Expected fewer number of vars in default output due to filtering but got $num_in and $num_out resp. (cmd was $cmd)"
fi


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


