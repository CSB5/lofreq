from distutils.core import setup, Extension
# see also http://docs.python.org/distutils/setupscript.html

import os
import sys
#import subprocess

import setup_conf


DEBUG = False
#DEBUG = True
                
DEFINE_MACROS = []
EXTRA_COMPILE_ARGS = []
if not DEBUG:
    DEFINE_MACROS.append(('NDEBUG', '1'))
    EXTRA_COMPILE_ARGS.append('-O3')
    EXTRA_COMPILE_ARGS.append('-funroll-loops')
else:
    EXTRA_COMPILE_ARGS.append('-Wall')
    EXTRA_COMPILE_ARGS.append('-Wextra')
    EXTRA_COMPILE_ARGS.append('-pedantic')
    EXTRA_COMPILE_ARGS.append('-ansi')
    EXTRA_COMPILE_ARGS.append('-O0')
    EXTRA_COMPILE_ARGS.append('-g')
    # NDEBUG seems to be always define through Python's OPT
    EXTRA_COMPILE_ARGS.append('-UNDEBUG') 


# C extensions
#
#PYAPI_PATH = os.path.join("src", "ext")
PYAPI_PATH = "lofreq_ext"
PYAPI_SOURCES = [os.path.join(PYAPI_PATH, f)
               for f in ["lofreq_ext.c"]]
CDF_LIBDIR = os.path.join('..', 'cdflib90')
CDF_LIB = os.path.join(CDF_LIBDIR, "libcdf.a")
LOFREQCORE_LIBDIR = os.path.join('..', 'lofreq_core')
LOFREQCORE_LIB = os.path.join(LOFREQCORE_LIBDIR, "liblofreq_core.a")


#def which(prog):
#    """make sure prog can be run
#    """
#    try:
#        subprocess.call([prog], 
#                        stderr=subprocess.PIPE, 
#                        stdout=subprocess.PIPE,)
#        return True
#    except OSError:
#        return False
#
    
# checks
#
if sys.version_info < (2 , 6):
    sys.stderr.write("FATAL: sorry, Python versions below 2.6 are not supported\n")
    sys.exit(1)
if sys.version_info >= (2 , 8):
    sys.stderr.write("FATAL: sorry, Python versions above 2.8 are not supported\n")
    sys.exit(1)

#for prog in []:
#    if not which(prog):
#        sys.stderr.write("#\nWARNING: cannot find '%s',"
#                         " which should have been installed earlier.\n#\n" % prog)
#        raw_input("Press Enter to continue anyway.")
#

extension = Extension("lofreq_ext",
                      PYAPI_SOURCES,
                      include_dirs=[CDF_LIBDIR, LOFREQCORE_LIBDIR],
                      define_macros=DEFINE_MACROS,
                      extra_compile_args=EXTRA_COMPILE_ARGS,
                      extra_objects=[CDF_LIB, LOFREQCORE_LIB],
                      depends=[CDF_LIB, LOFREQCORE_LIBDIR],
                      # libs statically linked using extra_objects instead of:
                      # libraries=['cdf'], library_dirs=[CDFLIBDIR],
                      )

# where modules reside:
package_dir = {'': setup_conf.PACKAGE_NAME.lower()}
    
setup(name = setup_conf.PACKAGE_NAME,
      packages=[setup_conf.PACKAGE_NAME.lower()],
      version = setup_conf.PACKAGE_VERSION,
      description="Low frequency variant caller",
      author="Andreas Wilm",
      author_email=setup_conf.PACKAGE_BUGREPORT,
      long_description = """LoFreq is a fast and sensitive variant-caller for inferring single-nucleotide variants (SNVs) from high-throughput sequencing data""",
      url='https://sourceforge.net/p/lofreq/',
      scripts = [
          'scripts/lofreq_alnoffset.py',
          'scripts/lofreq_bonf.py',
          'scripts/lofreq_diff.py',
          'scripts/lofreq_filter.py',
          'scripts/lofreq_pileup_summary.py',
          'scripts/lofreq_regionbed.py',
          'scripts/lofreq_snpcaller.py',
          'scripts/lofreq_uniq.py',
          'scripts/lofreq_uniq_pipeline.py',
          'scripts/lofreq_varpos_to_vcf.py',
          'scripts/lofreq_version.py'
      ],
      ext_modules = [extension],
      # http://pypi.python.org/pypi?%3Aaction=list_classifiers
      classifiers=['Development Status :: 4 - Beta',
                   'Environment :: Console',
                   'Intended Audience :: Science/Research',
                   'License :: OSI Approved :: GNU General Public License (GPL)',
                   'Natural Language :: English',
                   'Operating System :: Unix',
                   'Programming Language :: C',
                   'Programming Language :: Python :: 2.6',
                   'Programming Language :: Python :: 2.7',
                   'Topic :: Scientific/Engineering :: Bio-Informatics',
                   ],
      keywords='bioinformatics'
      )

