#!/bin/bash

source lib.sh || exit 1

basedir=../../lofreq-test-data//alnoffset/


# offsetting snp positions to self via alignment should give identical
# results

gotmd5=$(lofreq_alnoffset.py -s 'consensus-fake' -o - \
    -i $basedir/fake_paired_sample.snp -a $basedir/fake_aln.fa | \
    lofreq_alnoffset.py -a $basedir/fake_aln.fa -m 'consensus-fake'  -i - -o - | \
    cut -f 1-3 -d ' ' | $md5)
expmd5=$(cat $basedir/fake_paired_sample.snp | cut -f 1-3 -d ' ' | $md5)
if [ $gotmd5 != $expmd5 ]; then
    echoerror "SNV offsetting to self via alignment does not return original"
else
    echook "SNV offsetting to self via alignment"
fi

# converting back to the src positions (before cloning) should give identical positions
title="Converting back to pre-cloning positions"
#cat<<EOF
if lofreq_alnoffset.py -s 'consensus-fake' -m 'consensus-real' \
    -i $basedir/fake_paired_sample.snp -o - -a $basedir/fake_aln.fa \
    | tr ':' ' ' | awk '{if ($2!=$NF) {exit 1}}'; then
    echook $title
else
    echoerror $title
fi
#EOF


