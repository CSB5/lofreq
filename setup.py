from distutils.core import setup, Extension
import os

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
    EXTRA_COMPILE_ARGS.append('-Wall')
    EXTRA_COMPILE_ARGS.append('-Wextra')
    EXTRA_COMPILE_ARGS.append('-pedantic')
    EXTRA_COMPILE_ARGS.append('-ansi')
    EXTRA_COMPILE_ARGS.append('-O0')
    EXTRA_COMPILE_ARGS.append('-g')
    # NDEBUG seems to be always define through Python's OPT
    EXTRA_COMPILE_ARGS.append('-UNDEBUG') 


import sys

if sys.version_info < (2 , 6):
    sys.stderr.write("FATAL: sorry, Python versions below 2.6 are not supported\n")
    sys.exit(1)
if sys.version_info >= (2 , 8):
    sys.stderr.write("FATAL: sorry, Python versions above 2.8 are not supported\n")
    sys.exit(1)

    
# C extensions
#
#EXT_PATH = os.path.join("src", "ext")
EXT_PATH = "lofreq_ext"
EXT_SOURCES = [os.path.join(EXT_PATH, f)
               for f in ["lofreq_ext.c", "fet.c",
                         os.path.join("cdflib90", "dcdflib.c"),
                         os.path.join("cdflib90", "ipmpar.c"),
                         ]
               ]

extension = Extension("lofreq_ext",
                      EXT_SOURCES,
                      include_dirs=[os.path.join(EXT_PATH, 'cdflib90')],
                      define_macros=DEFINE_MACROS,
                      extra_compile_args=EXTRA_COMPILE_ARGS,
                      depends=[os.path.join(EXT_PATH, "cdflib90", "cdflib.h")],
                      #library_dirs=[],
                      #libraries=[],
                      )


# where modules reside:
package_dir = {'': 'lofreq'}
    
setup(name = "lofreq_snpcaller",
      packages=['lofreq'],
      version = "0.1",
      description="Low frequency variant caller",
      author="Andreas Wilm",
      author_email='wilma@gis.a-star.edu.sg',
      #long_description = """FIXME.""" 
      url='http://www.gis.a-star.edu.sg/', # FIXME
      scripts = ['scripts/lofreq_snpcaller.py',
                 #'scripts/lofreq_bonf.py',
                 'scripts/lofreq_pileup_summary.py',
                 #'scripts/lofreq_strandbias.py',
                 'scripts/lofreq_filter.py',
                 'scripts/lofreq_unique.py',
                 'scripts/lofreq_varpos_to_vcf.py'
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

