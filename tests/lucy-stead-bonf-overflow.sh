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

ref=/projects/delrosariorc/Bcells/Mapping_Reference/human_g1k_v37.fasta
bam=../../lofreq-test-data/lucy-stead-2012-11-08_smalltesting.bam
# used to crash with
# Traceback (most recent call last):
# ...
#   File "/home/ufaserv1_m/medlste/bin/LoFreq-0.3.2/lib/python2.7/site-packages/lofreq/qual.py", line 157, in call_snp_in_column
#     self.bonf_factor, self.sig_thresh)
# TypeError: Failed to parse arguments.
# because of Bonferroni factor overflow
lofreq_snpcaller.py -b $bam -f $ref --format vcf -o /dev/null --force || exit 1
echook "Run completed succesfully"
