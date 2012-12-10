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

# md5sum is md5 on mac
md5=$(which md5sum 2>/dev/null || which md5)

bam=../../lofreq-test-data/denv2-multiplex-replicates/ATCACG_1.bam
reffa=../../lofreq-test-data/denv2-multiplex-replicates/ref.fa
bed=../../lofreq-test-data/denv2-multiplex-replicates/region.bed

snv_out_bauto=${bam%.bam}_bauto.snp
snv_out_bset=${bam%.bam}_bset.snp
snv_out_b1=${bam%.bam}_b1.snp

# delete output files from previous run
DEBUG=0
if [ $DEBUG -ne 1 ]; then
	rm -f $snv_out_bset $snv_out_bauto $snv_out_b1 2>/dev/null
fi
# index bam if necessary
test -s ${bam}.bai || samtools index $bam


# determine bonferroni factor
lofreq_snpcaller.py --bonf $(lofreq_bonf.py --bam $bam) \
   -f $reffa -b $bam -o $snv_out_bset || exit 1
lofreq_snpcaller.py \
   -f $reffa -b $bam -o $snv_out_bauto || exit 1
lofreq_snpcaller.py --bonf 1\
   -f $reffa -b $bam -o $snv_out_b1 || exit 1

md5bset=$(cat $snv_out_bset | $md5)
md5bauto=$(cat $snv_out_bauto | $md5)
md5b1=$(cat $snv_out_b1 | $md5)

if [ "$md5bset" != "$md5bauto" ]; then
    echoerror "SNVs predicted using internally computed Bonf factor differ from the ones predicted with manually set factor"
elif [ "$md5bauto" == "$snv_out_b1" ]; then
    echoerror "SNVs predicted using internally computed Bonf factor should be different from using bonf 1"
else
    echook "Automatic Bonferroni settings work as expected"
fi
