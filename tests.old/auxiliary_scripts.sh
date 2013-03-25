#!/bin/bash

source lib.sh || exit 1


reffa=../../lofreq-test-data/denv2-multiplex-replicates/ref.fa
regionarg="-r consensus:1-100"
bam=../../lofreq-test-data/denv2-multiplex-replicates/ATCACG_1.bam
if lofreq_samtools mpileup -f $reffa $regionarg $bam 2>/dev/null | \
    lofreq_pileup_summary.py | tail -n 1 | grep -q ^100; then
    echook "lofreq_pileup_summary.py with lofreq_samtools works"
else
    echoerror "lofreq_pileup_summary.py with lofreq_samtools input failed"
fi
if samtools mpileup -f $reffa $regionarg $bam 2>/dev/null | \
    lofreq_pileup_summary.py --orig-samtools | tail -n 1 | grep -q ^100; then
    echook "lofreq_pileup_summary.py with samtools input works"
else
    echoerror "lofreq_pileup_summary.py with samtools input failed"
fi
    

bam=../../lofreq-test-data/denv2-simulation/denv2-10haplo.bam
reffa=../../lofreq-test-data/denv2-simulation/denv2-refseq.fa
num_var_pos=$(lofreq_varpos_to_vcf.py -t .05 -b $bam -o - -r $reffa | wc -l)
if [ $num_var_pos -lt 40 ]; then
    echoerror "Was more expecting variant positions determined by lofreq_varpos_to_vcf.py"
else
    echook "lofreq_varpos_to_vcf.py gave expected result"
fi
