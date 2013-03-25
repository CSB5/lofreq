#!/bin/bash

source lib.sh || exit 1

basedir=../../lofreq-test-data/denv2-ngc-diff-ref/
msa=$basedir/ref_aln.fa

refseq_bam=$basedir/GA004-SR-R00075_CTTGTA_s_2_denv2-refseq/vipr_bwa-uniq-dups-marked-recal.bam
refseq_diffsnp=$basedir/GA004-SR-R00075_CTTGTA_s_2_denv2-refseq/vipr_bwa-uniq-dups-marked-recal_lofreq-final_lowfreqonly_alnoffset_diff.snp
refseq_reffa=$basedir/den2_refseq.fa

eu660415_bam=$basedir/GA004-SR-R00075_CTTGTA_s_2_denv2-EU660415/vipr_bwa-uniq-dups-marked-recal.bam
eu660415_diffsnp=$basedir/GA004-SR-R00075_CTTGTA_s_2_denv2-EU660415/vipr_bwa-uniq-dups-marked-recal_lofreq-final_lowfreqonly_alnoffset_diff.snp
eu660415_reffa=$basedir/den2_EU660415.fa 

nuniq=$(lofreq_uniq.py -d $eu660415_diffsnp -b $refseq_bam \
    -r $refseq_reffa --chrom NC_001474.2 -a $msa -m NC_001474.2 | wc -l)
nexp=0
if [ $nuniq != $nexp ]; then
    echoerror "Expected $nexp but got $nuniq uniq SNVs"
else
    echook "Got expected number of uniq SNVs"
fi

# pos 40 to 59 likely undetected primer contamination
# no idea why there are four left
nuniq=$(lofreq_uniq.py -d $refseq_diffsnp \
    -b $eu660415_bam \
    -r $eu660415_reffa --chrom EU660415 -a $msa -m EU660415 | \
    awk '{if ($2<40 || $2>59) {print}}' | wc -l)
nexp=4
if [ $nuniq -gt $nexp ]; then
    echoerror "Expected at max $nexp but got $nuniq uniq SNVs"
else
    echook "Got expected number of uniq SNVs"
fi

# another way to test whether the positions are correct would be
#  to run in debug and check that ref always > snp at testing snps
