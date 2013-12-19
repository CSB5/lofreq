#!/bin/bash

source lib.sh || exit 1

bam=./data/denv2-pseudoclonal/denv2-pseudoclonal.bam
reffa=./data/denv2-pseudoclonal/denv2-pseudoclonal_cons.fa
bed=./data/denv2-pseudoclonal/denv2-pseudoclonal_incl.bed
vcf=$(mktemp -t $(basename $0).XXXXXX.vcf)
rm -f $vcf

# index bam if necessary
test -s ${bam}.bai || samtools index $bam

$LOFREQ call -f $reffa -l $bed -b auto -o $vcf $bam || exit 1
# this tests 'filter' as well as it's part of call
#export  PERL5LIB=/Users/wilma/local/lib/
#if perl -mVcf -e validate ../tests/denv2-pseudoclonal.vcf; then
if vcf-validator ../tests/denv2-pseudoclonal.vcf; then
    echook "Got valid VCF output"
else
    echoerror "Invalid VCF output"
fi
