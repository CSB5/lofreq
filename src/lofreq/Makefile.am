AM_CFLAGS=-D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE -Wall \
	-I../libbam/ -I../lofreq_core/ -I../uthash @AM_CFLAGS@
bin_PROGRAMS = lofreq
lofreq_SOURCES = defaults.h lofreq_bamstats.c lofrreq_bamstats.h lofreq_filter.c lofreq_filter.h lofreq_main.c lofreq_snpcaller.c lofreq_snpcaller.h lofreq_uniq.c lofreq_uniq.h lofreq_vcfset.h lofreq_vcfset.c
# note: order matters
lofreq_LDADD = ../lofreq_core/liblofreq_core.a ../libbam/libbam.a ../cdflib90/libcdf.a

