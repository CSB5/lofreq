#!/bin/bash

# FIXME:add-doc

source lib.sh || exit 1

basedir=data/denv2-simulation
bam=$basedir/denv2-10haplo.bam
reffa=$basedir/denv2-refseq.fa
truesnv=$basedir/denv2-10haplo_true-snp.snp

outdir=$(mktemp -d)
outraw=$outdir/raw.vcf
outfinal=$outdir/final.vcf

$LOFREQ2 call -b 30000 -f $reffa -o $outraw $bam || exit 1

$LOFREQ2 filter -p --strandbias-holmbonf --min-cov 10 -i $outraw -o $outfinal || exit 1


snvs_raw=$(cut -f1,2,4,5 $outraw | grep -v '^#' | tr '\t' ' ' | sort)
snvs_true=$(cut -f 1-3 -d ' ' $truesnv | tr '>' ' ' | sort

nexp=$(echo $snvs_raw | wc -l)
nraw=$(echo $snvs_true | wc -l)
if [ $nexp -ne $ncommon ]; then
    echoerror "Number of common SNVs differs (expected $nexp got $ncommon)"
fi

md5exp=$(echo $snvs_raw | $md5)
md5raw=$(echo $snvs_true | $md5)
if [ $md5exp != $md5raw ]; then
    echoerror "Number of SNV matches but content differs"
fi

echook "Tests passed"



# FIXME outfinal should not look different, i.e. filtering shouldn't do much/anything.
# see /home/wilma/snpcaller/lofreq/lofreq-sourceforge.git/tests/denv2-simulation.sh 

echoerror "Compare $outraw and $outfinal with $truesnv"
rm $outraw $outfinal
rmdir $outdir

