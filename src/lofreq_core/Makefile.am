AM_CFLAGS = -Wall -I../cdflib90/ -I../libbam/ -fPIC
AM_LDFLAGS = -shared -fPIC
noinst_LIBRARIES = liblofreq_core.a
liblofreq_core_a_SOURCES = bam2depth.c binom.c fet.c snpcaller.c utils.c
#liblofreq_core_a_LIBADD = ../cdflib90/libcdf.a ../libbam/libbam.a


#install-exec-hook:
#	 cd ../lofreq_python && LDFLAGS='-L$(prefix)/lib' $(PYTHON) setup.py install --prefix '$(prefix)'
#
#clean-local:
#	 cd ../lofreq_python && $(PYTHON) setup.py clean
#	 cd ../lofreq_python && rm -rf ./build

