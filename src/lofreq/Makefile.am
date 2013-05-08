AM_CFLAGS=-D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE -Wall \
	-I../libbam/ -I../lofreq_core/ @AM_CFLAGS@
bin_PROGRAMS = lofreq
lofreq_SOURCES = lofreq_main.c lofreq_uniq.c lofreq_uniq.h lofreq_snpcaller.c lofreq_snpcaller.h
# note: order matters
lofreq_LDADD = ../lofreq_core/liblofreq_core.a ../libbam/libbam.a ../cdflib90/libcdf.a

