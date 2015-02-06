#!/bin/bash

# running valgrind on call incl indel calls on ecoli spikein

source lib.sh || exit 1

KEEP_TMP=0
BASEDIR=data/ecoli-clone/
BAM=$BASEDIR/spike-in/spike-in.bam
REF=$BASEDIR/ref/Ecoli_K12_MG1655_NC_000913.fa

for f in $BAM $REF; do
  if [ ! -e $f ]; then
      echoerror "Required file $f missing"
      exit 1
  fi
done
            

outdir=$(mktemp -d -t $(basename $0).XXXXXX)
log=$outdir/log.txt
valgrindlog=$outdir/valgrind.log
vcf_out=$outdir/out.vcf


# how to get a region with true SNVs and indels close-by
#ipython
#import vcf
#vcfr = vcf.Reader(filename="truth.vcf.gz")
#vars = [v for v in vcfr]
#indel_highq = [v for v in vars_highq if v.is_indel]
#snv_highq = [v for v in vars_highq if v.is_snp]
#def argmin(iterable):
#        return min(enumerate(iterable), key=lambda x: x[1])[0]
#def closest(v, cmp_list):
#    dists = [abs(v.POS-c.POS) for c in cmp_list]
#    return argmin(dists)
#for i in indel_highq:
#    c = closest(i, snv_highq)
#    print i, snv_highq[c]
# and check that both are present in truth and lofreq prediction

valgrind --suppressions=faidx_fetch_seq.supp --leak-check=full --tool=memcheck --log-file=$valgrindlog \
  $LOFREQ call --call-indels -f $REF $BAM -r 'NC_000913:2000-2600' -o $vcf_out >$log 2>&1 || exit 1

for pos in 2000 2032 2214 2514 2572; do 
  if ! grep -q -w $pos $vcf_out; then
    echoerror "Excepted variant position $pos not found in vcf $vcf_out"
    exit 1
  fi
done
echook "All expected variant positions found"


num_err=$(grep 'ERROR SUMMARY' $valgrindlog | grep -cv ': 0 errors')
if [ "$num_err" -ne 0 ]; then
    echoerror "Found errors in Valgrind output $valgrindlog"
    exit 1
else
    echook "No errors found in Valgrind output"
fi

lost_bytes=$(grep 'lost' $valgrindlog | grep -cv ': 0 bytes in 0 blocks')
if [ "$lost_bytes" -ne 0 ]; then
    echoerror "Found lost bytes in Valgrind output $valgrindlog" || exit 1
    exit 1
else
    echook "No lost bytes found in Valgrind output"
fi

if [ $KEEP_TMP -ne 1 ] && [ $num_err -eq 0 ]; then
    test -d $outdir && rm -rf $outdir
else
    echowarn "Not deleting temporary output directory $outdir"
fi
