"""setup routines for LoFreq"""

from distutils.core import setup, Extension
from distutils.command.install import install as DistutilsInstall
import os
import sys
import subprocess

# see also http://docs.python.org/distutils/setupscript.html

DEBUG = False
#DEBUG = True
             
DEFINE_MACROS = []
EXTRA_COMPILE_ARGS = []
if not DEBUG:
    DEFINE_MACROS.append(('NDEBUG', '1'))
    EXTRA_COMPILE_ARGS.append('-O3')
    EXTRA_COMPILE_ARGS.append('-funroll-loops')
else:
    print "Note: Also consider setting DISTUTILS_DEBUG to anything (except the empty string)"
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
#EXT_PATH = os.path.join("src", "ext")
EXT_PATH = "lofreq_ext"
EXT_SOURCES = [os.path.join(EXT_PATH, f)
               for f in ["lofreq_ext.c", 
                         "fet.c",
                         os.path.join("cdflib90", "dcdflib.c"),
                         os.path.join("cdflib90", "ipmpar.c"),
                         ]
               ]


def which(prog):
    """make sure prog can be run
    """
        
    try:
        subprocess.call([prog], 
                        stderr=subprocess.PIPE, 
                        stdout=subprocess.PIPE,)
        return True
    except OSError:
        return False


def build_lofreq_samtools():
    """guess what: builds lofreq_samtools 
    """

    log_filename = 'Makefile.log'
    origdir = os.getcwd()
    
    try:
        os.chdir('lofreq_samtools')
        log = open(log_filename, 'w')
        subprocess.check_call(['make'], stdout=log, stderr=log)
        log.close()
    except:
        sys.stderr.write(
            "\nERROR: couldn't build lofreq_samtools!"
            " For more info, have a look at %s\n\n" % (
                os.path.join(os.getcwd(), log_filename)))
        sys.exit(1)
        
    os.chdir(origdir)
        
        
# custom install class so that we can install external source
# See http://stackoverflow.com/questions/1754966/how-can-i-run-a-makefile-in-setup-py
# 
class MyInstall(DistutilsInstall):
    def run(self):
        # pre install stuff goes here
        build_lofreq_samtools()

        DistutilsInstall.run(self)
       
        # post install stuff goes her
        if not which('lofreq_samtools'):
            sys.stderr.write("\nERROR: Looks like"
                " lofreq_samtools installation failed.")
            sys.exit(1)

        
sys.stderr.write("FIXME: test if lofreq_samtools is packed for dist as well\n")


# checks
#
if sys.version_info < (2 , 6):
    sys.stderr.write("FATAL: sorry, Python versions below 2.6 are not supported\n")
    sys.exit(1)
if sys.version_info >= (2 , 8):
    sys.stderr.write("FATAL: sorry, Python versions above 2.8 are not supported\n")
    sys.exit(1)

#for prog in ['samtools']:
#    if not which(prog):
#        sys.stderr.write("#\nWARNING: cannot find '%s',"
#                         " which is required to run LoFreq\n#\n" % prog)
#        raw_input("Press Enter to continue anyway.")


extension = Extension("lofreq_ext",
                      EXT_SOURCES,
                      include_dirs=[os.path.join(EXT_PATH, 'cdflib90')],
                      define_macros=DEFINE_MACROS,
                      extra_compile_args=EXTRA_COMPILE_ARGS,
                      depends=[os.path.join(EXT_PATH, "cdflib90", "cdflib.h"), 
                          os.path.join(EXT_PATH, "fet.h")]
                      #library_dirs=[],
                      #libraries=[],
                      )

# where modules reside:
package_dir = {'': 'lofreq'}
    
setup(name = 'LoFreq',
      packages=['lofreq'],
      version = '0.4.0',
      description="Low frequency variant caller",
      author="Andreas Wilm",
      author_email='wilma@gis.a-star.edu.sg',
      #long_description = """FIXME.""" 
      url='https://sourceforge.net/p/lofreq/', # FIXME
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
          'scripts/lofreq_version.py',
           # EVIL HACK! 'lofreq_samtools/lofreq_samtools'
      ],
      data_files=[
          # Also an evil hack. Wouldn't be surprised if this breaks
          # FIXME: test on Linux with and without prefix
          ('bin', ['lofreq_samtools/lofreq_samtools'])],
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
      keywords='bioinformatics',
      cmdclass={'install': MyInstall}
      )

