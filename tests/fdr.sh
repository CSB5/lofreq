#!/bin/bash

source lib.sh || exit 1

VCF=data/vcf/fdr.vcf
# 25 simulated variants extracted from mq-demo.
# converted probabilities given in fdr example in 
# http://www.biostathandbook.com/multiplecomparisons.html
# (see also multtest.c) to qualities and replace original values.


# expecting following result which mimicks the same as in link and in multtest.c
# expecting 5 significant results
NEXP=5
nres=$($LOFREQ filter --no-defaults -q fdr -r 0.25 -i $VCF | grep -vc '^#')
if [ $nres -ne $NEXP ]; then
    echoerror "FDR filtering not producing expected results (got $nres instead of $NEXP)"
    exit 1
fi

# even after capping and setting #tests
vcf=$(mktemp -t $(basename $0).XXXXXX.vcf)
head -n 11 $VCF > $vcf
nres=$($LOFREQ filter --no-defaults -q fdr -r 0.25 -s 25 -i $vcf | grep -vc '^#')
if [ $nres -ne $NEXP ]; then
    echoerror "FDR filtering after capping not producing expected results (got $nres instead of $NEXP)"
    exit 1
fi
echook "FDR filtering produced expected results"
