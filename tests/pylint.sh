#!/bin/bash

scripts2test=$(grep 'scripts/' ../setup.py | tr -d ",'" | sed -e 's,^,../,')
log=$(mktemp -t pylint.XXXXX)
for f in $scripts2test; do
    pylint -E --rcfile pylint.rc $f >> $log
done 
if [ -s $log ]; then
    cat $log
    exit 1
else
    echo "pylint produced no errors"
fi
rm $log


