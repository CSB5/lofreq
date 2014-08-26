#!/bin/bash

source lib.sh || exit 1

vcf=data/icgc-tcga-dream-testproject/strelka-1.0.13_snvs-indels-somatic.vcf

num_total=$(grep -vc '^#' $vcf)
num_snvs=$($LOFREQ filter --no-default --only-snvs -i $vcf | grep -vc '^#')
num_indels=$($LOFREQ filter --no-default --only-indels -i $vcf | grep -vc '^#')
vcf=tests/data/icgc-tcga-dream-testproject/strelka-1.0.13_snvs-indels-somatic.vcf
msg="Number of SNVs ($num_snvs) and indels ($num_indels) extracted by filter"
if [ $(expr $num_snvs + $num_indels) -ne $num_total ]; then
     echoerror "$msg don't add up to total number of variants ($num_total)"
     exit 1
else
     echook "$msg add up to total number of variants ($num_total)"
fi
