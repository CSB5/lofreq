# Top-level Makefile.am for LoFreq
#
# based on Makefile.am for Clustal Omega

# we do not use GNU-style files: NEWS README AUTHORS ChangeLog
# AUTOMAKE_OPTIONS = foreign

SUBDIRS = src/lofreq_samtools src/lofreq_python

bug-tests: all
	cd tests && $(SHELL) run_all.sh;

