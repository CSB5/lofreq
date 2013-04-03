#!/bin/bash

# FIXME:add-doc

source lib.sh || exit 1


vcf_t=data/vcf/CHH966-tumor-100x-100pur-hg19.bwa_6431925.vcf.gz
vcf_n=data/vcf/CHH966-normal-100x-100pur-hg19.bwa.renamed_6431925.vcf.gz
vcf_out=/tmp/FIXMEtest.vcf

../src/lofreq_python/scripts/lofreq2_vcfset.py \
  -1 $vcf_t -2 $vcf_n -a complement -o - | cut -f 1-7 > $vcf_out

num_diffs=$(gzip -dc data/vcf/CHH966-tumor-only.f-7.vcf.gz | \
    diff -u $vcf_out - | grep -v '##' | grep '^[\+\-]' | wc -l)
exp_diffs=10
#--- test.vcf	2013-04-03 22:12:53.000000000 +0800
#+++ -	2013-04-03 22:22:06.000000000 +0800
#-chr12	30805918	.	C	G	23	.
#-chr13	107516488	.	T	G	22	.
#-chr16	69170707	.	G	C	23	.
#-chr17	8738690	.	T	G	23	.
#-chr2	42513376	.	G	C	23	.
#-chr4	186560162	.	C	G	22	.
#-chr6	42571331	.	T	A	24	.
#-chr6	106553829	.	G	A	26	.
#
# All diffs expected. vcf-isec only looks at chrom and pos, not the 

if [ $num_diffs -ne $exp_diffs ]; then
    echoerror "Expected $exp_diffs but got $num_diffs (keeping $vcf_out for your reference)."
else
    echook "Complement test produced expected results."
    rm $vcf_out
fi

