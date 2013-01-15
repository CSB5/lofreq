AM_CFLAGS=-D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE -Wall \
	-I../libbam/ -L../libbam/
bin_PROGRAMS = lofreq_samtools
lofreq_samtools_SOURCES = bam_plcmd.c
lofreq_samtools_LDADD = ../libbam/libbam.a

