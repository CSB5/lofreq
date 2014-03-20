#!/bin/bash

myname=$(basename $0)

source lib.sh || exit 1

PY_DIRS="../src/tools/lofreq_star"
files=$(find $PY_DIRS -name \*py -not -name _\*)
for f in $files; do
    echo "$myname: testing $f"
    python $f || echoerror "testing $f failed"
done

for f in $files; do
    echo "$myname: testing $f"
	python -m doctest $f || echoerror "testing $f failed"
done
