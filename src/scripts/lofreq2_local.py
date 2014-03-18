# add local dir to path to make source dir, i.e. not installed scripts
# work straight-away

import sys
import os

# Set sys.path/PYTHONPATH such that we find the local source dir first
# by using: from lofreq_star import ...
#d = os.path.normpath(os.path.join(
#    os.path.dirname(sys.argv[0]), '..'))
#if os.path.exists(os.path.join(d, "lofreq_star")):
#    #sys.stderr.write("NOTE: Adding local dir %s to PYTHONPATH\n" % d)
#    sys.path.insert(0, d)

# Set PATH such that we find lofreq binary first
d = os.path.normpath(os.path.join(
    os.path.dirname(sys.argv[0]), '../lofreq'))
if os.path.exists(os.path.join(d, 'lofreq')):
    #sys.stderr.write("NOTE: Adding local dir %s to PATH\n" % d)
    os.environ["PATH"] = d + os.pathsep + os.environ["PATH"]

# In theory need to find scripts because the main binary knows about them. However, there are circular cases where script call the binary which then can't find the scripts again (e.g. in parallel wrapper),so:
#
#d = os.path.normpath(os.path.join(
#    os.path.dirname(sys.argv[0]), '../tools/scripts'))
#if os.path.exists(d):
#    #sys.stderr.write("NOTE: Adding local dir %s to PATH\n" % d)
#    os.environ["PATH"] = d + os.pathsep + os.environ["PATH"]
    

        
