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

# index bam if necessary
test -s ${bam}.bai || samtools index $bam

echowarn "Only checking that something's output but should check validity"
lofreq_snpcaller.py -f $reffa -l $bed -b $bam -o /dev/null --force --format vcf || exit 1
echook "Got some VCF output"
