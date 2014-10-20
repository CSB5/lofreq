# Top-level Makefile.am for LoFreq
#
# based on Makefile.am for Clustal Omega

# we do not use GNU-style files: NEWS README AUTHORS ChangeLog
# AUTOMAKE_OPTIONS = foreign

if ENABLE_TOOLS
TOOLS = src/tools
endif

LICENSES = LICENSE src/cdflib90.README src/zlib.LICENSE src/libbam/AUTHORS src/libbam/COPYING src/uthash/LICENSE
EXTRA_DIST = binary_installer.sh $(LICENSES)

SUBDIRS = src/cdflib90 src/libbam src/lofreq src/scripts src/uthash src/lofreq_alnqual $(TOOLS) 
bug-tests: all
	cd tests && $(SHELL) run_all.sh;

