#!/bin/bash

source lib.sh || exit 1

KEEP_TMP=0
REF=./data/af_tests/ref_fasta.fa
outdir=$(mktemp -d -t $(basename $0).XXXXXX)

# See ./data/af_tests/README for expected results
failed=0

# del test
bam=./data/af_tests/test_deletions.bam
log=$outdir/del_log.txt
vcf=$outdir/del_out.vcf
cmd="$LOFREQ call --no-default-filter -B -f $REF -o $vcf $bam"
#echodebug "cmd=$cmd"
if ! eval $cmd > $log 2>&1; then
    echoerror "LoFreq failed. Check logfile $log. Command was $cmd"
    exit 1
fi
if ! awk '{if ($2=="1" && $4=="ACG" && $5=="A" && $8 ~ /AF=0.18/) {m=1; exit 0}} END {if (m) {exit 0} else {exit 1}}' $vcf; then
    echoerror "Expected deletion of AF=0.5 not found in $vcf"
    let failed=failed+1
fi
if ! awk '{if ($2=="1" && $4=="A" && $5=="T" && $8 ~ /AF=1.0/) {m=1; exit 0}} END {if (m) {exit 0} else {exit 1}}' $vcf; then
    echoerror "Expected SNV of AF=1.0 not found in $vcf"
    let failed=failed+1
fi

# ins test
bam=./data/af_tests/test_insertion.bam
log=$outdir/ins_log.txt
vcf=$outdir/ins_out.vcf
cmd="$LOFREQ call --no-default-filter -B -a 0.5 -f $REF -o $vcf $bam"
#echodebug "cmd=$cmd"
if ! eval $cmd > $log 2>&1; then
    echoerror "LoFreq failed. Check logfile $log. Command was $cmd"
    exit 1
fi
if ! awk '{if ($2=="2" && $4=="C" && $5=="CAA" && $8 ~ /AF=0.5/) {m=1; exit 0}} END {if (m) {exit 0} else {exit 1}}' $vcf; then
    echoerror "Expected insertion of AF=1.0 not found in $vcf"
    let failed=failed+1
fi
if ! awk '{if ($2=="2" && $4=="C" && $5=="G" && $8 ~ /AF=0.25/) {m=1; exit 0}} END {if (m) {exit 0} else {exit 1}}' $vcf; then
    echoerror "Expected SNV of AF=0.25 not found in $vcf"
    let failed=failed+1
fi


# FIXME check output
if [ $KEEP_TMP -ne 1 ] && [ $failed -eq 0 ]; then
   test -d $outdir && rm -rf $outdir
fi
                    
                    

