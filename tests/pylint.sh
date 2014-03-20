#!/bin/bash

source lib.sh || exit 1

#files_to_test=$(grep 'scripts/' ../src/lofreq_python/setup.py | tr -d "[,']" | tr -d '[\t ]' | sed -e 's,^,../,')
#files_to_test=$(grep '^[^#].*\.py'  ../src/lofreq_python/Makefile.am | grep -v PYTHON | cut -f 2 -d = | tr -d '\\' | tr -d '[\t ]')
files_to_test=$(find ../src/scripts ../src/tools/scripts ../src/tools/lofreq_star -name \*py)
PYLINT=$(which pylint 2>/dev/null || which pylint-2.7) || exit 1

echoinfo "Using $PYLINT"
log=$(mktemp -t pylint.XXXXX)
for f in $files_to_test; do
    echoinfo "Testing $f"
    $PYLINT -E --rcfile pylint.rc $f >> $log
done 
if [ -s $log ]; then
    echoerror "pylint produced errors:"
    cat $log
    exit 1
else
    echook "pylint produced no errors"
fi
rm $log


