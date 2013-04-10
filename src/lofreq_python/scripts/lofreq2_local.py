# Add local dir to path to make source dir, i.e. not installed scripts
# work straight-away

import sys, os

d = os.path.join(
    os.path.dirname(sys.argv[0]), '..')
if os.path.exists(os.path.join(d, 'lofreq2')):
    #sys.stderr.write("Adding local dir %s to PYTHONPATH\n" % d)
    sys.path.insert(0, d)

        
