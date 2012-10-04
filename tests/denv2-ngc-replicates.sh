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


# in
bam_1=../example-data/denv2-ngc-replicates/GGCTAC_2_remap_razers-i92_peakrem_corr.bam
ref_1=../example-data/denv2-ngc-replicates/GGCTAC_2_remap_razers-i92_peakrem.covcons.fa
snp_1=../example-data/denv2-ngc-replicates/GGCTAC_2_lofreq.snp
# created here
snp_sbf_1=../example-data/denv2-ngc-replicates/GGCTAC_2_lofreq_sbf.snp
snp_diff_1=../example-data/denv2-ngc-replicates/GGCTAC_2_lofreq_sbf_diff-raw.snp

# in
bam_2=../example-data/denv2-ngc-replicates/CTTGTA_2_remap_razers-i92_peakrem_corr.bam 
ref_2=../example-data/denv2-ngc-replicates/CTTGTA_2_remap_razers-i92_peakrem.covcons.fa
snp_2=../example-data/denv2-ngc-replicates/CTTGTA_2_lofreq.snp
# create here
snp_sbf_2=../example-data/denv2-ngc-replicates/CTTGTA_2_lofreq_sbf.snp
snp_diff_2=../example-data/denv2-ngc-replicates/CTTGTA_2_lofreq_sbf_diff-raw.snp


# delete output files from previous run
rm -f $snp_sbf_1 $snp_diff_1 $snp_sbf_2 $snp_diff_2 2>/dev/null

# index bam if necessary
for bam in $bam_1 $bam_2; do
    test -s ${bam}.bai || samtools index $bam
done



# sb-filtering
#
echo -e "$snp_1 $snp_sbf_1 58 1st\n$snp_2 $snp_sbf_2 38 2nd" | \
    while read fsnp fsnp_sbf nexp what; do
    lofreq_filter.py --strandbias-holmbonf -i $fsnp > $fsnp_sbf
    nsbfsnvs=$(wc -l $fsnp_sbf | cut -f 1 -d ' ')
    
    if [ $nexp -ne $nsbfsnvs ]; then
        echoerror "Number of sb-filtered SNVs differs (expected $nexp got $nsbfsnvs) for $what file"
    else
        echook "Got expected number of sb-filtered SNVs for $what file"
    fi
done


# diffing
#
lofreq_diff.py -s $snp_sbf_1 -t $snp_sbf_2 -m uniq_to_1 > $snp_diff_1
ndiff=$(wc -l $snp_diff_1 | cut -f 1 -d ' ')
nexp=30
if [ $nexp -ne $ndiff ]; then
    echoerror "Number of raw diff SNVs in 1st file differ"
else
    echook "Got expected number of raw diff SNVs in 1st file"
fi

lofreq_diff.py -s $snp_sbf_1 -t $snp_sbf_2 -m uniq_to_2 > $snp_diff_2
ndiff=$(wc -l $snp_diff_2 | cut -f 1 -d ' ')
nexp=10
if [ $nexp -ne $ndiff ]; then
    echoerror "Number of raw diff SNVs in 2nd file differ"
else
    echook "Got expected number of raw diff SNVs in 2nd file"
fi


# unique sanity (against themselves)
#
echo -e "$snp_diff_1 $ref_1 $bam_1 1st\n$snp_diff_2 $ref_2 $bam_2 2nd" | \
    while read diff ref bam what; do
    nuniq=$(lofreq_unique.py -d $diff -r $ref -b $bam | grep -c ': unique')
    if [ 0 -ne $nuniq ]; then
        echoerror "Number of unique SNVs in sanity check on $what file differ (got $nuniq instead of 0)"
    else
        echook "Got expected number (0) unique SNVs in sanity check on $what file"
    fi
done



# unique
#
echo -e "$snp_diff_1 $ref_2 $bam_2 7 1st\n$snp_diff_2 $ref_1 $bam_1 8 2nd" | \
    while read diff ref bam nexp what; do
    nuniq=$(lofreq_unique.py -d $diff -r $ref -b $bam | grep -c ': unique')
    if [ $nexp -ne $nuniq ]; then
        echoerror "Number of unique SNVs differ first $what diff (got $nuniq expected $nexp)"
    else
        echook "Got expected number of unique SNVs in $what diff"
    fi
done


