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


bam1=../example-data/denv2-ngc-replicates/CTTGTA_2_remap_razers-i92_peakrem_corr.bam
bam2=../example-data/denv2-ngc-replicates/GGCTAC_2_remap_razers-i92_peakrem_corr.bam
# choose either as ref (they are identical)
ref=../example-data/denv2-ngc-replicates/GGCTAC_2_remap_razers-i92_peakrem.covcons.fa
outdir=../example-data/denv2-ngc-replicates-uniq-pipeline/

bed=$outdir/region.bed
# keep original snv calls but delete all others
rm -f $(ls $outdir/*snp 2>/dev/null | grep -v raw.snp) >/dev/null
rm $bed 2>/dev/null

lofreq_regionbed.py -i $bam1 > $bed || exit 1
lofreq_uniq_pipeline.py --bam1 $bam1 --bam2 $bam2 \
    --bed $bed --ref $ref -o $outdir || exit 1

# lazy: just check total number of uniq snvs
nuniq=$(cat $outdir/*uniq.snp | wc -l)
nexp=11
if [ $nexp -ne $nuniq ]; then
    echoerror "Expected number of uniq SNV differs" 
else
    echook "Got expected number of uniq SNVs"
fi


