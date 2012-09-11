find lofreq/ -name \*py -not -name _\* | \
	xargs -n1 python
find scripts/ -name \*py -not -name _\* | \
	xargs -n1 python -m doctest 


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
