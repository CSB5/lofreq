#!/bin/bash

# FIXME:add-doc

source lib.sh || exit 1


vcf_t=data/vcf/CHH966-tumor-100x-100pur-hg19.bwa_6431925.vcf.gz
vcf_n=data/vcf/CHH966-normal-100x-100pur-hg19.bwa.renamed_6431925.vcf.gz
#vcf_t=data/vcf/CHH966-tumor-100x-100pur-hg19.bwa_6431925.vcf
#vcf_n=data/vcf/CHH966-normal-100x-100pur-hg19.bwa.renamed_6431925.vcf
vcf_out=$(mktemp -t $(basename $0).XXXXXX.vcf)

cmd="$LOFREQ vcfset -1 $vcf_t -2 $vcf_n -a complement -o -"

#echodebug "cmd=$cmd"
eval $cmd | cut -f 1-7 > $vcf_out

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
    echook "Larger complement test produced expected results."
    rm $vcf_out
fi





vcf_1=data/vcf/vcf_set.vcf.gz
vcf_1_allfiltered=data/vcf/vcf_set_allfiltered.vcf.gz

# complement against self should give zero
cmd="$LOFREQ vcfset -1 $vcf_1 -2 $vcf_1 -a complement -o -"
num_compl=$(eval $cmd | grep -vc '^#')
if [ $num_compl -ne 0 ]; then
    echoerror "Complement against self should give 0"
else
    echook "Complement against self returned 0"
fi


# intersect against self should give all
cmd="$LOFREQ vcfset -1 $vcf_1 -2 $vcf_1 -a intersect -o -"
md5_test=$(eval $cmd | grep -v '^#' | $md5)
md5_org=$(zgrep -v '^#' $vcf_1 | $md5)
if [ "$md5_test" != "$md5_org" ]; then
    echoerror "Intersect against self should give results identical to input (cmd: $cmd)"
    #echodebug "md5_test = $md5_test"
    #echodebug "md5_org = $md5_org"
else
    echook "Intersect against self gave results identical to input"
fi


# intersect with all filtered should give 0
cmd="$LOFREQ vcfset -1 $vcf_1 -2 $vcf_1_allfiltered -a intersect --only-passed -o -"
num_inter=$(eval $cmd | grep -vc '^#')
if [ $num_inter -ne 0 ]; then
    echoerror "Intersect with all filtered should give 0 (but gave $num_inter; cmd = $cmd)"
else
    echook "only-passed Intersect with all filtered returned 0"
fi

# complement with all filtered should give all
cmd="$LOFREQ vcfset -1 $vcf_1 -2 $vcf_1_allfiltered --only-passed -a complement -o -"
md5_test=$(eval $cmd | grep -v '^#' | $md5)
md5_org=$(zgrep -v '^#' $vcf_1 | grep 'PASS' | $md5)
#echodebug "$cmd test=$md5_test org=$md5_org"
if [ "$md5_test" != "$md5_org" ]; then
    echoerror "only-passed complement with all filtered should give results identical to input (cmd = $cmd)"
else
    echook "only-passed complement with all filtered gave results identical to input"
fi


#
vcf_org=data/vcf/vcf_set.vcf.gz
vcf_baseswap=data/vcf/vcf_set_altrefswap.vcf.gz
cmd="$LOFREQ vcfset -1 $vcf_org -2 $vcf_baseswap -a intersect -o -"
num_out=$(eval $cmd | grep -cv '^#')
if [ $num_out -ne 0 ]; then
    echoerror "intersection with base swapped file did not return any variants"    
else
    echook "intersection with base swapped file return variants"    
fi

cmd="$cmd --only-pos"
num_out=$(eval $cmd | grep -cv '^#')
if [ $num_out -eq 0 ]; then
    echoerror "intersection with base swapped file when using bases did not return zero variants"
else
    echook "intersection with base swapped file return zero variants"    
fi

