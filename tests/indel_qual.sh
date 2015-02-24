#!/bin/bash

source lib.sh || exit 1

set -o pipefail

# test correct use of indel qualities

sam=data/idq/delq_3lq.sam
ref=data/idq/ref.fa
res=$(samtools view -bS $sam 2>/dev/null | $LOFREQ call --call-indels -f $ref  -b 1 -a 1 --no-default-filter -B -A - 2>&1 || exit 1)
#echo "$res"
if ! echo "$res" | grep -q 'ref[[:space:]]1'; then
	echoerror "Should have called indel at pos 1 but didn't (res was $res)"
	exit 1
fi
if echo "$res" | grep -q 'ref[[:space:]]3'; then
        echoerror "Shouldn't have called indel at pos 3 but did (res was $res)"
        exit 1
fi

sam=data/idq/delq_1lq.sam
res=$(samtools view -bS $sam 2>/dev/null | $LOFREQ call --call-indels -f $ref  -b 1 -a 1 --no-default-filter -B -A - 2>&1 || exit 1)
#echo "$res"
if echo "$res" | grep -q 'ref[[:space:]]1'; then
        echoerror "Shouldn't have called indel at pos 1 but did (res was $res)"
        exit 1
fi
if ! echo "$res" | grep -q 'ref[[:space:]]3'; then
       echoerror "Should have called indel at pos 3 but didn't (res was $res)"
       exit 1
fi

echook "Indels predicted as expected"