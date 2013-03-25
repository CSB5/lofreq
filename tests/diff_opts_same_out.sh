#!/bin/bash

# Test to make sure that paralell computed results are identical to
# the ones computed without paralllel option and that reading from
# stdin also results in same result

source lib.sh || exit 1


echowarn "Better use a less high coverage data-set for faster completion"
indir=../../lofreq-test-data/denv2-pseudoclonal/
bam=$indir/denv2-pseudoclonal.bam
ref=$indir/denv2-pseudoclonal_cons.fa
bed=$indir/denv2-pseudoclonal_incl.fake.bed


CMD[1]="$LOFREQ2 call -l $bed -f $ref --verbose $bam"
CMD[2]="$LOFREQ2 call -l $bed -f $ref --verbose --pseudo-parallel 4 $bam"
CMD[3]="cat $bam | $LOFREQ2 call -f $ref -l $bed --verbose -"
# cannot: CMD[4]="cat $bam | $LOFREQ2 call -f $ref -l $bed --verbose --pseudo-parallel 4 -"
for i in $($seq 1 ${#CMD[@]}); do
    cmd=${CMD[$i]}
    out=$(mktemp -t $(basename $0).XXXXXX.vcf)
    log=$(mktemp -t $(basename $0).XXXXXX.log)

    #echodebug "Executing $cmd with output going to $out and $log"
    # remove source line from vcf which will change depending on call
    if ! eval $cmd 2>$log | grep -v 'source' >$out ; then
        echoerror "Executing following command failed (see $log for more info): $cmd"
        exit 1
    fi

    # make sure we predicted at least one snv. if output is always
    # empty tests would be successful otherwise
    if ! grep -q DP4 $out; then
        echoerror "No SNVs in output file $out found"
        exit 1
    fi

    # compare to output of previous cmd
    if [ -n "$prevout" ]; then
        if ! diff -q $out $prevout; then
            echoerror "Results between runs differed. Commands were:"
            echoerror " Current cmd:  $cmd"
            echoerror " Current out:  $out"
            echoerror " Previous cmd: $prevcmd"
            echoerror " Previous out: $prevout"
            exit 1
        fi
    fi
    
    prevcmd=$cmd
    test -s "$prevlog" && rm $prevlog
    test -s "$prevout" && rm $prevout
    prevlog=$log
    prevout=$out
done
test -s "$prevlog" && rm $prevlog
test -s "$prevout" && rm $prevout






