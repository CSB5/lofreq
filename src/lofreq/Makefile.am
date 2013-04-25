AM_CFLAGS=-O2 -D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE -Wall \
	-I../libbam/ -I../lofreq_core/
bin_PROGRAMS = lofreq
lofreq_SOURCES = lofreq_main.c lofreq_snpcaller.c lofreq_snpcaller.h
lofreq_LDADD = ../lofreq_core/liblofreq_core.a ../libbam/libbam.a

