#!/bin/bash

source lib.sh || exit 1

bam=../../lofreq-test-data/denv2-multiplex-replicates/ATCACG_1.bam
reffa=../../lofreq-test-data/denv2-multiplex-replicates/ref.fa
bed=../../lofreq-test-data/denv2-multiplex-replicates/region.bed

snv_out_bauto=${bam%.bam}_bauto.snp
snv_out_bautonon0=${bam%.bam}_bautonon0.snp
snv_out_bset=${bam%.bam}_bset.snp
snv_out_b1=${bam%.bam}_b1.snp

# delete output files from previous run
DEBUG=0
if [ $DEBUG -ne 1 ]; then
	rm -f $snv_out_bset $snv_out_bauto $snv_out_b1 $snv_out_bautonon0 2>/dev/null
fi
# index bam if necessary
test -s ${bam}.bai || samtools index $bam


# determine bonferroni factor
lofreq_snpcaller.py --bonf $(lofreq_bonf.py --bam $bam) \
   -f $reffa -b $bam -o $snv_out_bset || exit 1
lofreq_snpcaller.py \
   -f $reffa -b $bam -o $snv_out_bauto || exit 1
lofreq_snpcaller.py \
   -f $reffa -b $bam -o $snv_out_bautonon0 --bonf auto-ign-zero-cov || exit 1
lofreq_snpcaller.py --bonf 1\
   -f $reffa -b $bam -o $snv_out_b1 || exit 1

md5bset=$(cat $snv_out_bset | $md5)
md5bauto=$(cat $snv_out_bauto | $md5)
md5bautonon0=$(cat $snv_out_bautonon0 | $md5)
md5b1=$(cat $snv_out_b1 | $md5)

if [ "$md5bset" != "$md5bauto" ]; then
    echoerror "SNVs predicted using internally computed Bonf factor differ from the ones predicted with manually set factor"
elif  [ "$md5bauto" != "$md5bautonon0" ]; then
    # only two zero depth pos
    echoerror "SNVs predicted using internally computed Bonf factor differ from the ones predicted with auto non-zero-cov factor"    
elif [ "$md5bauto" == "$snv_out_b1" ]; then
    echoerror "SNVs predicted using internally computed Bonf factor should be different from using bonf 1"
else
    echook "Automatic Bonferroni settings work as expected"
fi
