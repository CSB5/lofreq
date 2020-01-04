# configure.ac for LoFreq

AC_PREREQ(2.63)
# 2.60 needed for AC_USE_SYSTEM_EXTENSIONS
# 2.64 which allows to define a URL as well
# 2.68 seems to have updated ax_pthread

AC_INIT([LoFreq_Star], [2.1.4], [wilma@gis.a-star.edu.sg])

# The AC_INIT macro can take any source file as an argument. It just
# checks that the file is there, which should, in turn, mean that the
# source directory is there. AC_INIT generates the following m4
# macros, output variables and preprocessor symbols:
AC_SUBST(PACKAGE_NAME)
AC_SUBST(PACKAGE_TARNAME)
AC_SUBST(PACKAGE_VERSION)
AC_SUBST(PACKAGE_STRING)
AC_SUBST(PACKAGE_BUGREPORT)

#AC_REVISION([m4_esyscmd_s([git describe --always --dirty])])
AC_REVISION([m4_esyscmd_s([git describe --always])])
# saves revision to configure

git_version=`git describe --always --dirty`
if test x"$git_version" = x; then
   AC_DEFINE(GIT_VERSION, ["unknown"])
else
   AC_DEFINE_UNQUOTED(GIT_VERSION, "$git_version")
fi

AC_MSG_NOTICE([Configuring $PACKAGE_NAME (version $PACKAGE_VERSION) for your system...

            |             ____|                 
            |       _ \   |     __|  _ \   _` | 
            |      (   |  __|  |     __/  (   | 
           _____| \___/  _|   _|   \___| \__, | 
                                             _| 
])
#
# Created with http://patorjk.com/software/taag/#p=display&f=Shadow&t=LoFreq


# The AM_INIT_AUTOMAKE line adds several standard checks. This macro
# is always needed for automake. Obsolete: It takes the program name
# and version number as arguments.
AM_INIT_AUTOMAKE([foreign])


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
# AC_PROG_RANLIB rendered obsolte by LT_INIT?
AC_PROG_RANLIB

#
# disable-static: Need to disable static libraries, otherwise we can't
# link them into the dynamically loaded Python extension
#
#LT_INIT([disable-static])

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
# Note, optimization level can be overriden by user with env var
# CFLAGS
AC_PROG_CC
# Set output variable CPP to a command that runs the C preprocessor. 
AC_PROG_CPP

# If the C compiler cannot compile ISO Standard C (currently C99),
# try to add an option to output variable `CC' to make it work.  If
# the compiler does not support C99, fall back to supporting ANSI
# C89 (ISO C90).
# LoFreq: In theory we need for variadic macro support 
# AC_PROG_CC_STDC
#
#If the C compiler is not in C99 mode by default, try to add an
#option to output variable `CC' to make it so.  This macro tries
#various options that select C99 on some system or another,
#preferring extended functionality modes over strict conformance
#modes.  It considers the compiler to be in C99 mode if it handles
#`_Bool', `//' comments, flexible array members, `inline', signed
#and unsigned `long long int', mixed code and declarations, named
#initialization of structs, `restrict', `va_copy', varargs macros,
#variable declarations in `for' loops, and variable length arrays.
#After calling this macro you can check whether the C compiler has
#been set to accept C99; if not, the shell variable
#`ac_cv_prog_cc_c99' is set to `no'.
AC_PROG_CC_C99



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

# pthread support
#
AC_MSG_CHECKING([if we have pthreads support])
m4_include([m4/ax_pthread.m4])
AC_CONFIG_MACRO_DIR([m4])
AX_PTHREAD([
        AC_DEFINE(HAVE_PTHREAD,1,[Define if you have POSIX threads libraries and header files.])
        CLIBS="$PTHREAD_LIBS $LIBS"
        CFLAGS="$CFLAGS $PTHREAD_CFLAGS"
        LDFLAGS="$LDFLAGS $PTHREAD_CFLAGS"
        CC="$PTHREAD_CC"],
        AC_MSG_ERROR([No pthread support on this machine]))
#AX_PTHREAD()


# explicit libm check
AC_CHECK_LIB(m, log,, AC_MSG_ERROR([Could not find libm]))
AC_CHECK_LIB(z, gzread,, AC_MSG_ERROR([Could not find libz]))


# http://www.gnu.org/software/automake/manual/html_node/Python.html
AM_PATH_PYTHON([2.7])
# argparse and a lot of other stuff is not available in 2.7

# on by default
AC_ARG_ENABLE([tools],
    AS_HELP_STRING([--disable-tools], [Disable installation of tools]))
AM_CONDITIONAL(ENABLE_TOOLS, [test "x$enable_tools" != "xno"])

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
    CFLAGS="${CFLAGS} -O0 -pedantic"
    CXXFLAGS="${CXXFLAGS} -O0 -pedantic"
else
    # switch off assert etc
    AC_DEFINE([NDEBUG],[],[No-debug Mode])
    # Could set -O3 but that should be a user choice (env var CFLAGS)
fi

if test x"$HTSLIB" != x""; then
   AC_MSG_ERROR([Use --with-htslib=/fullpath/to/htslibdir instead of HTSLIB variable])
fi
m4_include([m4/ax_with_htslib.m4])
AX_WITH_HTSLIB
if test "$ax_cv_htslib" != yes; then
   AC_MSG_ERROR([HTSlib not found])
fi
AM_CONDITIONAL([HTSLIB_IS_SRC], [test "$ax_cv_htslib_which" = source])

AC_SUBST([AM_CFLAGS])
AC_SUBST([AM_LDFLAGS])

# AC_DEFINE([SOURCEQUAL_IGNORES_INDELS], [1], [ignore indels in SQ computation as long as we can't predict them])


#AC_ARG_ENABLE([alnerrprof],
#        [AS_HELP_STRING([--enable-alnerrprof], [Enable use of alignment/mapping error profile (default=no)])],
#        AC_DEFINE([USE_ALNERRPROF], [1], [support mapping error profile in LoFreq]))


AC_CONFIG_FILES(Makefile
        src/cdflib90/Makefile   
        src/uthash/Makefile
        src/lofreq/Makefile
        src/scripts/Makefile        
        src/tools/setup_conf.py        
        src/tools/Makefile
        )
AC_OUTPUT

AC_MSG_NOTICE([Done. Next, run 'make' to compile the source, followed by 'make install' to install it.])
