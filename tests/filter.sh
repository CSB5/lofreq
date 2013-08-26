#!/bin/bash

# Test that we get the number of expected SNVs on the pseudo-clonal data-set

source lib.sh || exit 1


vcf=data/vcf/CTTGTA_2_remap_razers-i92_peakrem_corr_nodeff.vcf
#outvcf=$(mktemp -t $(basename $0).XXXXXX)

base_cmd="$LOFREQ filter -i $vcf --no-defaults -o -"

ALPHA_LIST='0.01 0.0001 0.000001 0.00000001'
NUMTEST_LIST='100 10000 1000000'


# snv quality with varying alpha
#
for cor in "bonf" "holm-bonf" "fdr"; do
    last_no=0
    for a in $ALPHA_LIST; do
        cmd="$base_cmd --snv-qual $cor --snv-qual-alpha $a"
        #echodebug "cmd=$cmd"
        new_no=$(eval $cmd | grep -c 'snvqual.*\(bonf\|fdr\)') || exit 1
        #echodebug "$cor a=$a: new_no=$new_no last_no=$last_no";# cmd = $cmd"
        if [ $new_no -lt $last_no ]; then
            echoerror "snvqual: Got fewer SNVs filtered with higher alpha (cmd=$cmd)"
        fi
        last_no=$new_no
    done
done

# snv quality with varying num_tests
#
# fixed alpha
a=0.00000001
for cor in "bonf" "holm-bonf" "fdr"; do
    last_no=0
    for n in $NUMTEST_LIST; do
        cmd="$base_cmd --snv-qual $cor --snv-qual-alpha $a --snv-qual-numtests $n"
        #echodebug "cmd=$cmd"
        new_no=$(eval $cmd | grep -c 'snvqual.*\(bonf\|fdr\)') || exit 1
        #echodebug "$cor a=$a n=$n: new_no=$new_no last_no=$last_no";# cmd = $cmd"
        if [ $new_no -lt $last_no ]; then
            echoerror "snvqual: Got fewer SNVs filtered with higher num-tests (cmd=$cmd)"
        fi
        last_no=$new_no
    done
done


# strandbias quality with varying alpha
#
for cor in "bonf" "holm-bonf"; do
    last_no=100000
    for a in $ALPHA_LIST; do
        cmd="$base_cmd --strandbias $cor --strandbias-alpha $a"
        #echodebug "cmd=$cmd"
        new_no=$(eval $cmd | grep -c 'strandbias.*bonf') || exit 1
        #echodebug "$cor a=$a: new_no=$new_no last_no=$last_no";# cmd = $cmd"
        if [ $new_no -gt $last_no ]; then
            echoerror "strandbias: Got more SNVs filtered with higher alpha (cmd=$cmd)"
        fi
        last_no=$new_no
    done
done
