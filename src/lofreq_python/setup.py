from distutils.core import setup
# see also http://docs.python.org/distutils/setupscript.html

import os
import sys
#import subprocess

import setup_conf

DEBUG = False
#DEBUG = True


# checks
#
if sys.version_info < (2 , 6):
    sys.stderr.write("FATAL: sorry, Python versions"
                     " below 2.6 are not supported\n")
    sys.exit(1)
if sys.version_info >= (2 , 8):
    sys.stderr.write("FATAL: sorry, Python versions"
                     " above 2.8 are not supported\n")
    sys.exit(1)

   
# where modules reside:
#package_dir = {'': setup_conf.PACKAGE_NAME.lower()}
#package_dir = {'': ''}

    
setup(name = setup_conf.PACKAGE_NAME,
      packages=[setup_conf.PACKAGE_NAME.lower()],
      version = setup_conf.PACKAGE_VERSION,
      description="Low frequency variant caller",
      author="Andreas Wilm",
      author_email=setup_conf.PACKAGE_BUGREPORT,
      long_description = """LoFreq-Star is a fast and sensitive variant-caller for inferring single-nucleotide variants (SNVs) from high-throughput sequencing data""",
      #url='https://sourceforge.net/p/lofreq/',
      scripts = [
          'scripts/lofreq2_filter.py',
          'scripts/lofreq2_somatic.py', 
          'scripts/lofreq2_uniq.py',
          'scripts/lofreq2_vcfset.py'
      ],
      # http://pypi.python.org/pypi?%3Aaction=list_classifiers
      classifiers=['Environment :: Console',
                   'Intended Audience :: Science/Research',
#'License :: OSI Approved :: GNU General Public License (GPL)',
                   'Natural Language :: English',
                   'Operating System :: Unix',
                   'Programming Language :: C',
                   'Programming Language :: Python :: 2.7',
                   'Topic :: Scientific/Engineering :: Bio-Informatics',
                   ],
      keywords='bioinformatics'
      )
