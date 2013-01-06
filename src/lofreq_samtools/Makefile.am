AM_CFLAGS=-D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE -D_USE_KNETFILE -Wall \
	-I../libbam/ -L../libbam/
bin_PROGRAMS = lofreq_samtools
lofreq_samtools_SOURCES = bam_plcmd.c bam2depth.c sam_view.c bamtk.c
lofreq_samtools_LDADD = ../libbam/libbam.la

