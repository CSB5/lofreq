# Add local dir to path to make source dir, i.e. not installed scripts
# work straight-away

import sys
import os

d = os.path.normpath(os.path.join(
    os.path.dirname(sys.argv[0]), '..'))
if os.path.exists(os.path.join(d, 'lofreq_star')):
    sys.stderr.write("NOTE: Adding local dir %s to PYTHONPATH\n" % d)
    sys.path.insert(0, d)

d = os.path.normpath(os.path.join(
    os.path.dirname(sys.argv[0]), '../../lofreq'))
if os.path.exists(os.path.join(d, 'lofreq')):
    sys.stderr.write("NOTE: Adding local dir %s to PATH\n" % d)
    os.environ["PATH"] = os.pathsep + d + os.environ["PATH"]
    

        
