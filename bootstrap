#!/bin/sh

# A bootstrapping script replacing autogen.sh and autoreconf. Brings
# source tree into a state where end user can run configure and make.
#
# From https://www.sourceware.org/autobook/autobook/autobook_43.html:
# "Autoconf comes with a program called autoreconf which essentially
# does the work of the bootstrap script. autoreconf is rarely used
# because, historically, has not been very well known, and only in
# Autoconf 2.13 did it acquire the ability to work with Automake.
# Unfortunately, even the Autoconf 2.13 autoreconf does not handle
# libtoolize and some automake-related options that are frequently
# nice to use.
#
# We recommend the bootstrap method, until autoreconf is fixed. At
# this point bootstrap has not been standardized, so here is a version
# of the script we used while writing this book"
#

aclocal && \
	automake --gnu --add-missing && \
	autoconf

