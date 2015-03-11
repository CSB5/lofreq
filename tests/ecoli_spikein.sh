#!/bin/bash

source lib.sh || exit 1

KEEP_TMP=1
BASEDIR=data/ecoli-clone/
BAM=$BASEDIR/spike-in/spike-in.bam
REF=$BASEDIR/ref/Ecoli_K12_MG1655_NC_000913.fa
TRUTH=$BASEDIR/spike-in/truth.laln.vcf.gz
EVALUATOR=data/icgc-tcga-dream-support/evaluator.py

for f in $BAM $REF $TRUTH $EVALUATOR; do
  if [ ! -e $f ]; then
    echoerror "Required file $f missing"
    exit 1
  fi
done

 
outdir=$(mktemp -d -t $(basename $0).XXXXXX)
log=$outdir/log.txt
outvcf=$outdir/out.vcf


cmd="$LOFREQ call-parallel --pp-threads 8 -f $REF -o $outvcf --verbose $BAM"
# only needed as long as indels are disabled by default
cmd="$cmd --call-indels"
echodebug "cmd=$cmd"
if ! eval $cmd > $log 2>&1; then
    echoerror "LoFreq failed. Check logfile $log. Command was $cmd"
    exit 1
fi


res_ll=$($EVALUATOR -v $outvcf -t $TRUTH -m SNV | awk 'END {print $NF}') || exit 1
res=$(echo $res_ll | \
  awk -F, '{prec=$1; rec=$2; if (prec<0.96 || rec<0.96) {status="ERROR"} else {status="OK"} printf "%s: precision=%f recall=%f\n", status, prec, rec}') || exit 1
if echo $res | grep -q ERROR; then
   let num_err=num_err+1
fi
echo "snvs: $res" 1>&2


# based on results for v2.1.2a-69-g5bd5919 -1%
res_ll=$($EVALUATOR -f $REF -v $outvcf -t $TRUTH -m INDEL | awk 'END {print $NF}') || exit 1
res=$(echo $res_ll | \
  awk -F, '{prec=$1; rec=$2; if (prec<0.956 || rec<0.917) {status="ERROR"} else {status="OK"} printf "%s: precision=%f recall=%f\n", status, prec, rec}') || exit 1
if echo $res | grep -q ERROR; then
   let num_err=num_err+1
fi
echo "indels: $res" 1>&2



if [ $KEEP_TMP -ne 1 ] && [ $num_err -eq 0 ]; then
    test -d $outdir && rm -rf $outdir
else
    echowarn "Not deleting temporary output directory $outdir"
fi
if [ $num_err -ne 0 ]; then
    exit 1
fi
