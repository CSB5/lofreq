#!/bin/bash

set -o pipefail

PREFIX=/usr/local

usage() {
    # keep in sync with arg parsing below
cat <<EOF
$(basename $0): binary installer for LoFreq

  Options:
    -p | --prefix  : Install to this directory (default: $PREFIX)
    -h | --help    : Display this help
EOF
}

while [ "$1" != "" ]; do
    case $1 in
        -p | --prefix )
            shift
            prefix=$1
            ;;
        -h | --help )
            usage
            exit
            ;;
        * ) 
            echo "FATAL: unknown argument \"$1\""
            usage
            exit 1
    esac
    shift
done

test -z "$prefix" && prefix=$PREFIX
echo "Using $prefix as installation prefix"

echo "Installing binaries"
test -d "$prefix/bin" || mkdir -p $prefix/bin || exit 1
BINARIES="./src/lofreq/lofreq src/scripts/lofreq2_call_pparallel.py src/scripts/lofreq2_somatic.py src/lofreq_alnqual/lofreq2_alnqual"
for f in $BINARIES; do
    if [ ! -s $f ]; then
        echo "FATAL: can't find $f" 1>&2
        exit 1
    fi
    cp $f $prefix/bin/ || exit 1
done

echo "Installing Python tools"
pushd ./src/tools >/dev/null
pythonpath=$(python setup.py install --prefix $prefix | \
    grep 'egg-info$' | head -n1 | cut -f 2 -d ' ' | sed -e 's,LoFreq.*,,') || exit 1

echo "NOTE: Make sure $pythonpath is in your PYTHONPATH"
popd >/dev/null

echo "Successful exit"
