#!/bin/bash

# Test that we get the number of expected SNVs on the pseudo-clonal data-set

source lib.sh || exit 1


VCF=data/vcf/filter_test.vcf.gz

FILTER="$LOFREQ filter --sb-no-compound"
#FILTER=../src/lofreq/lofreq_filter

# must be ordered
ALPHA_LIST='0.01 0.00001 0.00000001'
MTC_TYPES='bonf holmbonf fdr'
NUMTEST_LIST='10000 1000000'

# number of input variants
num_in=$(zgrep -vc '^#' $VCF)

# number of failed tests
num_fails=0


# #input == #output with and without filtering
#
num_out=$($FILTER -i $VCF -v 1 -V 2 -a 0.5 -A 0.6 -B 10 -q bonf | grep -vc '^#')
if [ $num_in -ne $num_out ]; then
    echoerror "total #input != #output (with filter)"
    let num_fails=num_fails+1
fi
num_out=$($FILTER -i $VCF | grep -vc '^#')
if [ $num_in -ne $num_out ]; then
    echoerror "total #input != #output (without filter)"
    let num_fails=num_fails+1
fi


# check defaults
#
num_filter_tags=$($FILTER -i $VCF | grep -v '^#' | grep -v PASS | cut -f 7 | tr ';' '\n' | sort -u | wc -l)
if [ $num_filter_tags -ne 2 ]; then
    echoerror "was expecting exactly two filter tags coming from default filtering"
    let num_fails=num_fails+1
fi


# AF filtering
#
num_exp=$(zgrep -c 'AF=0.2' $VCF)
num_out=$($FILTER -i $VCF --no-defaults --af-min 0.2 --af-max 0.3 --only-passed | grep -vc '^#')
if [ $num_exp -ne $num_out ]; then
    echoerror "AF filtering failed"
    let num_fails=num_fails+1
fi


# DP filtering
#
num_exp=$(zgrep -c 'DP=2[0-9][0-9]' $VCF)
num_out=$($FILTER -i $VCF --no-defaults  --cov-min 200 --cov-max 300 --only-passed | grep -vc '^#')
if [ $num_exp -ne $num_out ]; then
    echoerror "DP filtering failed"
    let num_fails=num_fails+1
fi


# SB threshold filtering
#
num_exp=$(zgrep -c 'SB=[0-9]\($\|;\)' $VCF)
num_out=$($FILTER -i $VCF --no-defaults  --sb-thresh 9 --only-passed | grep -vc '^#')
if [ $num_exp -ne $num_out ]; then
    echoerror "SB thresholdfiltering failed"
    let num_fails=num_fails+1
fi


# SB MTC
#
num_prev_mtc=100000
for mtc in $MTC_TYPES; do
    num_prev_alpha=$($FILTER -i $VCF --sb-mtc $mtc --only-passed | grep -vc '^#')

    # bonf rejects fewer than holm-bonf than fdr, i.e. in that order
    # more are significant, i.e. are filtered and therefore fewer pass
    #
    if [ $num_prev_alpha -gt $num_prev_mtc ]; then
        echoerror "SB $mtc produced let more variants pass then previous one"
        let num_fails=num_fails+1
        break
    fi
    num_prev_mtc=$num_prev_alpha
    
    # as alpha goes up, we become more stringent, i.e. fewer are
    # significant and more pass
    #
    for alpha in $ALPHA_LIST; do
        num_higher_alpha=$($FILTER -i $VCF --sb-mtc $mtc --sb-alpha $alpha --only-passed | grep -vc '^#')
        #echodebug "$mtc  $alpha  $num_prev_alpha -> $num_higher_alpha"
        if [ $num_higher_alpha -lt $num_prev_alpha ]; then
            echoerror "SB $mtc with next highest alpha ($alpha) produced fewer PASSED variants"
            let num_fails=num_fails+1
            break
        fi
        num_prev_alpha=$num_higher_alpha

    done
done



# SNV qual threshold filtering
#
Q=40
num_exp=$($zcat $VCF | awk -v q=$Q '/^[^#]/ {if ($6=="." || $6>=q) {s+=1}} END {print s}')
num_out=$($FILTER -i $VCF --no-defaults  --snvqual-thresh $Q --only-passed | grep -vc '^#')
if [ $num_exp -ne $num_out ]; then
    echoerror "SNV quality threshold filtering failed: expected $num_exp but got $num_out"
    let num_fails=num_fails+1
fi

if [ $num_fails -gt 0 ];then
    echoerror "$num_fails tests failed"
else
    echook "all tests passed"
fi


# SB MTC
#
num_prev_mtc=0
for mtc in $MTC_TYPES; do
    num_prev_alpha=$($FILTER -i $VCF --no-defaults --snvqual-mtc $mtc --only-passed | grep -vc '^#')

    # bonf rejects fewer than holm-bonf than fdr, i.e. in that order
    # more are significant, i.e. are kept and therefore more pass
    #
    if [ $num_prev_alpha -lt $num_prev_mtc ]; then
        echoerror "SNV qual $mtc produced let fewer variants pass than previous one"
        let num_fails=num_fails+1
        break
    fi
    num_prev_mtc=$num_prev_alpha
    
    # as alpha goes up, we become more stringent, i.e. fewer
    # significant and fewer pass
    #
    for alpha in $ALPHA_LIST; do
        num_higher_alpha=$($FILTER -i $VCF --no-defaults --snvqual-mtc $mtc --snvqual-alpha $alpha --only-passed | grep -vc '^#')
        #echodebug "$mtc  $alpha  $num_prev_alpha -> $num_higher_alpha"
        if [ $num_higher_alpha -gt $num_prev_alpha ]; then
            echoerror "SNV qual $mtc with next highest alpha ($alpha) produced more PASSED variants"
            let num_fails=num_fails+1
            break
        fi
        num_prev_alpha=$num_higher_alpha

    done
done


echo "WARN: manual diff against py impl. missing" 1>&2
# see:
# VCF_IN=data/vcf/filter_test.vcf.gz
# ../src/lofreq/lofreq filter -i $VCF_IN \
#     --no-defaults --cov-min 10 --cov-max \
#     --af-min 0.1 \
#     --sb-mtc bonf --sb-alpha 0.001\
#     --snvqual-mtc holmbonf --snvqual-alpha 1 --snvqual-ntests 100000 | \
#     grep -v '^#' | sed -e 's,_dp,cov,' -e 's,_,,g' -e 's,sb,strandbias,' | \
#     grep -c PASS
# 
# ../src/lofreq_python/scripts/lofreq2_filter.py -i $VCF_IN \
#     --no-defaults -o - \
#     --min-cov 10 --max-cov 90 \
#     --min-af 0.1 \
#     --strandbias bonf --strandbias-alpha 0.001 \
#     --snv-qual holm-bonf --snv-qual-alpha 1 --snv-qual-numtests 100000 | \
#     grep -v '^#' | \
#     grep -c PASS
# 
