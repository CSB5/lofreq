# configure.ac for LoFreq
#
# Based on configure.ac for Clustal Omega


AC_PREREQ(2.6) # 2.60 needed for AC_USE_SYSTEM_EXTENSIONS
# consider 2.64 which allows to define a URL as well

# FIXME set version here and replace in setup.py and samtools Makefile
# or any .h
AC_INIT([Lofreq], [0.5.0a], [wilma@gis.a-star.edu.sg])

# The AC_INIT macro can take any source file as an argument. It just
# checks that the file is there, which should, in turn, mean that the
# source directory is there. AC_INIT generates the following m4
# macros, output variables and preprocessor symbols:
AC_SUBST(PACKAGE_NAME)
AC_SUBST(PACKAGE_TARNAME)
AC_SUBST(PACKAGE_VERSION)
AC_SUBST(PACKAGE_STRING)
AC_SUBST(PACKAGE_BUGREPORT)

AC_MSG_NOTICE([configuring $PACKAGE_NAME (version $PACKAGE_VERSION) for your system...])


# The AM_INIT_AUTOMAKE line adds several standard checks. This macro
# is always needed for automake. Obsolete: It takes the program name
# and version number as arguments.
AM_INIT_AUTOMAKE

# config header: provides a way for people to customise the
# configuration for their platform, via #defines. config.h.in can be
# generated automatically with the autoheader tool. However, you need
# a stamp-h file in your project to ensure that automake regenerates
# config.h from config.h.in. Type 'touch stamp-h' to add this file to
# your project.
#
#AM_CONFIG_HEADER(src/config.h) 

# the following allows to install a copy config.h for libclustalo in a
## namespace like way, i.e. no name clash with your own config.h
# AX_PREFIX_CONFIG_H(src/clustal-omega-config.h)

# This macro was introduced in Autoconf 2.60. If possible, enable
# extensions to C or Posix on hosts that normally disable the
# extensions.
AC_USE_SYSTEM_EXTENSIONS


# for cross compilation
#
AC_CANONICAL_HOST
## AC_CANONICAL_BUILD
## AC_CANONICAL_TARGET

# Set the @SET_MAKE@ variable=make if $(MAKE) not set
AC_PROG_MAKE_SET

# AC_PROG_RANLIB will allow you to build code in sub-directories into
# temporary libraries, which make will then link in with the rest of
# the code.
# AC_PROG_RANLIB rendered obsolte by LT_INIT
#
# disable-static: Need to disable static libraries, otherwise we can't
# link them into the dynamically loaded Python extension
#
LT_INIT([disable-static])

# AC_PROG_INSTALL will generate an install target so that users may
# just type 'make install' to install the software.
AC_PROG_INSTALL

# check for C++ preprocessor and compiler and the library compiler
# (might change the compiler flags)
#
# Determine a C compiler to use. If using the GNU C compiler, set
# shell variable GCC to `yes'. If output variable CFLAGS was not
# already set, set it to `-g -O2' for the GNU C compiler (`-O2' on
# systems where GCC does not accept `-g'), or `-g' for other
# compilers.
AC_PROG_CC
# Set output variable CPP to a command that runs the C preprocessor. 
AC_PROG_CPP
#
# Determine a C++ compiler to use. If output variable CXXFLAGS was not
# already set, set it to `-g -O2' for the GNU C++ compiler (`-O2' on
# systems where G++ does not accept `-g'), or `-g' for other
# compilers.
#AC_PROG_CXX
# Set CXXCPP to a command that runs the C++ preprocessor
#AC_PROG_CXXCPP
#
# reset compiler flags to initial flags
#CXXFLAGS=$OLD_CXXFLAGS

# safety check: AC_CONFIG_SRCDIR([src/clustal.h])


# Checks for standard header files.
# obsolescent according to http://www.gnu.org/software/hello/manual/autoconf/Particular-Headers.html
AC_HEADER_STDC
# Guess the next one is obsolete as well then?
# FIXME AC_CHECK_HEADERS(assert.h ctype.h float.h limits.h math.h memory.h stdarg.h stdio.h stdlib.h string.h time.h unistd.h)
#  awk  '/^#include *</ {print $2}' $(find src/ -name \*.[ch] -or -name \*.cpp) | grep -v argtable2.h  | sort -u | tr -d '<>'

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_TYPE_SIZE_T
AC_C_INLINE

# Checks for library functions.
# AC_CHECK_FUNCS(strdup strstr strchr erand48)
# AW is there an automatic way to extract library functions used in the program?


# FIXME explicit libm check
#AC_CHECK_LIB(m, log,, AC_MSG_ERROR([Could not find libm]))
#AC_CHECK_LIB(m, log2, AC_DEFINE(HAVE_LOG2, 1, [Has log2()]))
AC_CHECK_LIB(m, log,, AC_MSG_ERROR([Could not find libm]))
AC_CHECK_LIB(z, gzread,, AC_MSG_ERROR([Could not find libz]))

# ---    argtable header and library check
#
#AC_CHECK_HEADERS([argtable2.h], [],
#                 AC_MSG_ERROR([Could not find argtable2.h. Try $ ./configure CFLAGS='-Iyour-argtable2-include-path]))
#AC_CHECK_LIB(argtable2, arg_parse, [],
#             AC_MSG_ERROR([Could not find libargtable2. Try $ ./configure LDFLAGS="-Lyour-argtable2-lib-path']))


# http://www.gnu.org/software/automake/manual/html_node/Python.html
AM_PATH_PYTHON([2.6])
#AC_MSG_CHECKING([for python-config to get compiler flags])
#AC_CHECK_PROG(python_config, 'python-config', yes)
#if test "${python_config}" != "yes"; then
#   AC_MSG_ERROR("Could not run python-config.")
#fi
#PYTHON_CFLAGS=`python-config --cflags`
#PYTHON_LDFLAGS=`python-config --ldflags`
#AC_SUBST(PYTHON_CFLAGS)
#AC_SUBST(PYTHON_LDFLAGS)



# ---   debugging
#
AC_MSG_CHECKING([whether to build with debug information])
AC_ARG_ENABLE([debug],
              [AS_HELP_STRING([--enable-debug],
                              [Enable debugging (default=no)])],
              [debugit="$enableval"],
              [debugit=no])
AC_MSG_RESULT([$debugit])
if test x"$debugit" = x"yes"; then
    # User variables CFLAGS/CXXFLAGS are used last during compilation and
    # can therefore overwrite system settings.
    #
    # if debugging was requested, add debugging info at the end of
    # C[XX]FLAGS, which makes sure previous optimisation flags are
    # overwritten (normally C[XX]FLAGS should not be modified as they
    # are user variables):
    # switch off assert etc
    # AC_DEFINE([DEBUG],[],[Debug Mode])
    CFLAGS="${CFLAGS} -g -ansi -Wall -O0"
    CXXFLAGS="${CXXFLAGS} -g -ansi -Wall -O0"
else
    # switch off assert etc
    AC_DEFINE([NDEBUG],[],[No-debug Mode])
    # if no debugging is requested, we make optimisation explicit
    AM_CXXFLAGS="${AM_CXXFLAGS} -O3"
    AM_CFLAGS="${AM_CFLAGS} -O3"
fi


AC_SUBST([AM_CFLAGS])
AC_SUBST([AM_LDFLAGS])

AC_DEFINE([_USE_KNETFILE],[],[KNETFILE for samtools])

AC_CONFIG_FILES(Makefile
        src/libbam/Makefile
        src/cdflib90/Makefile
        src/lofreq_core/Makefile
        src/lofreq_samtools/Makefile
        src/lofreq_python/setup_conf.py        
        src/lofreq_python/Makefile)
AC_OUTPUT
