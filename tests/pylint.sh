#!/bin/bash

pylint=$(which pylint 2>/dev/null || which pylint-2.7)

echoerror() {
    echo "ERROR: $1" 1>&2
}
echook() {
    echo "OK: $1" 1>&2
}

scripts2test=$(grep 'scripts/' ../setup.py | tr -d ",'" | sed -e 's,^,../,')
log=$(mktemp -t pylint.XXXXX)
for f in $scripts2test; do
    $pylint -E --rcfile pylint.rc $f >> $log
done 
if [ -s $log ]; then
    echoerror "pylint produced errors:"
    cat $log
    exit 1
else
    echook "pylint produced no errors"
fi
rm $log


