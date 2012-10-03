#!/bin/bash

myname=$(basename $0)

errorecho() {
    echo "ERROR: $@" 1>&2
}

for f in $(find ../lofreq/ -name \*py \
    -not -name _\* -not -name \*20[0-9][0-9]\*); do
    echo "$myname: testing $f"
    python $f || errorecho "testing $f failed"
done

for f in $(find ../scripts/ -name \*py \
    -not -name _\* -not -name \*20[0-9][0-9]\*); do
    echo "$myname: testing $f"
	python -m doctest $f || errorecho "testing $f failed"
done


#VERBOSE = False
#
#if __name__ == "__main__":
#    
#    for f in itertools.chain(
#        glob.glob('./lofreq/*py'),
#        glob.glob('./scripts/*py')
#        ):
#        if os.path.basename(f).startswith("_"):
#            continue
#        print "Testing %s" % f
#        doctest.testfile(f, VERBOSE)
#        print    
