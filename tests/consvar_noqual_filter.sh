#!/bin/bash

# FIXME:add-doc

source lib.sh || exit 1

VCF=data/vcf/consvar_only.vcf.gz

num_in=$(zgrep -vc '^#' $VCF)
num_out=$($LOFREQ filter --snvqual-thresh 1 --no-defaults --only-passed -i $VCF | grep -vc '^#')
if [ $num_in -ne $num_out ]; then
    echoerror "Some CONSVARs were filtered by snvqual-thresh."
else
    echook "CONSVARs untouched by snvqual-thresh filtering"
fi

