#!/bin/bash

source lib.sh || exit 1

bam=../../lofreq-test-data/denv2-pseudoclonal/denv2-pseudoclonal.bam
reffa=../../lofreq-test-data/denv2-pseudoclonal/denv2-pseudoclonal_cons.fa
bed=../../lofreq-test-data/denv2-pseudoclonal/denv2-pseudoclonal_incl.bed
vcf=$(basename $bam .bam).vcf

# index bam if necessary
test -s ${bam}.bai || samtools index $bam

echowarn "Only checking that something's output but should check validity"
lofreq_snpcaller.py -f $reffa -l $bed -b $bam -o $vcf --force --format vcf || exit 1
#export  PERL5LIB=/Users/wilma/local/lib/
if perl -mVcf -e validate ../tests/denv2-pseudoclonal.vcf; then
    echook "Got valid VCF output"
else
    echoerror "Invalid VCF output"
fi
