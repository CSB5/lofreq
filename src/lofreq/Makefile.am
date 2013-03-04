AM_CFLAGS=-D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE -Wall \
	-I../libbam/ -I../lofreq_core/
bin_PROGRAMS = lofreq2
lofreq2_SOURCES = lofreq_main.c lofreq_snpcaller.c
lofreq2_LDADD = ../libbam/libbam.a ../lofreq_core/liblofreq_core.a

