#!/bin/bash

# Test that we get the number of expected SNVs on the pseudo-clonal data-set

source lib.sh || exit 1


vcf=data/vcf/CTTGTA_2_remap_razers-i92_peakrem_corr_nodeff.vcf.gz
#outvcf=$(mktemp -t $(basename $0).XXXXXX)

# FIXME base_cmd="$LOFREQ filter -i $vcf --no-defaults -o -"
#base_cmd="../src/lofreq/lofreq_filter -i $vcf -o -"
base_cmd="$LOFREQ filter -i $vcf -o -"

ALPHA_LIST='0.01 0.0001 0.000001 0.00000001'
NUMTEST_LIST='100 10000 1000000'



# snv quality with varying alpha
#
num_fail=0
for cor in "bonf" "holm-bonf" "fdr"; do
    last_no=0
    for a in $ALPHA_LIST; do
        #cmd="$base_cmd --snv-qual $cor --snv-qual-alpha $a"
        cmd="$base_cmd --snvqual-mtc $cor --snvqual-alpha $a"
        #echodebug "cmd=$cmd"
        new_no=$(eval $cmd | grep -c 'snvqual.*\(bonf\|fdr\)') || exit 1
        #echodebug "$cor a=$a: new_no=$new_no last_no=$last_no";# cmd = $cmd"
        if [ $new_no -lt $last_no ]; then
            echoerror "snvqual: Got fewer SNVs when filtering with higher alpha (cmd=$cmd)"
            let num_fail=num_fail+1
        fi
        last_no=$new_no
    done
done
if [ $num_fail -eq 0 ]; then
    echook "snvqual (var alpha): all tests passed"
fi


# snv quality with varying num_tests
#
# fixed alpha
a=0.00000001
for cor in "bonf" "holm-bonf" "fdr"; do
    last_no=0
    for n in $NUMTEST_LIST; do
        #cmd="$base_cmd --snv-qual $cor --snv-qual-alpha $a --snv-qual-numtests $n"
        cmd="$base_cmd --snvqual-mtc $cor --snvqual-alpha $a --snvqual-ntests $n"
        #echodebug "cmd=$cmd"
        new_no=$(eval $cmd | grep -c 'snvqual.*\(bonf\|fdr\)') || exit 1
        #echodebug "$cor a=$a n=$n: new_no=$new_no last_no=$last_no";# cmd = $cmd"
        if [ $new_no -lt $last_no ]; then
            echoerror "snvqual: Got fewer SNVs when filtering with higher num-tests (cmd=$cmd)"
            let num_fail=num_fail+1
        fi
        last_no=$new_no
    done
done
if [ $num_fail -eq 0 ]; then
    echook "snvqual (var num_tests): all tests passed"
fi


# strandbias quality with varying alpha
#
num_fail=0
for cor in "bonf" "holm-bonf"; do
    last_no=100000
    for a in $ALPHA_LIST; do
        #cmd="$base_cmd --strandbias $cor --strandbias-alpha $a"
        cmd="$base_cmd --sb-mtc $cor --sb-alpha $a"
        #echodebug "cmd=$cmd"
        new_no=$(eval $cmd | grep -c 'strandbias.*bonf') || exit 1
        #echodebug "$cor a=$a: new_no=$new_no last_no=$last_no";# cmd = $cmd"
        if [ $new_no -gt $last_no ]; then
            echoerror "strandbias: Got more SNVs when filtering with higher alpha (cmd=$cmd)"
            let num_fail=num_fail+1
        fi
        last_no=$new_no
    done
done
if [ $num_fail -eq 0 ]; then
    echook "strandbias: all tests passed"
fi

# window filter
# FIXME: not implemented in C version
#
#num_fail=0
#base_cmd="$LOFREQ filter -i $vcf --no-defaults -o -"
#cmd="$base_cmd --window 10"
#num_reg=$(eval $cmd | grep '[^0-9,]85' | grep -c snvwin) || exit 1
#num_exp=4
#if [ $num_reg -ne $num_exp ]; then
#    echoerror "window: Got $num_reg but expected $num_exp SNVs (cmd = $cmd)"
#    let num_fail=num_fail+1
#fi
##
#cmd="$base_cmd --window 1"
#num_reg=$(eval $cmd | grep '[^0-9,]85' | grep -c snvwin) || exit 1
#num_exp=2
#if [ $num_reg -ne $num_exp ]; then
#    echoerror "window: Got $num_reg but expected $num_exp SNVs (cmd = $cmd)"
#    let num_fail=num_fail+1
#fi
#if [ $num_fail -eq 0 ]; then
#    echook "window: all tests passed"
#fi
#
