# Add these guys as well during 'make dist'
# find . -name \*py -or -name \*.c -or -name \*.h | grep -v /build/ | grep -v setup_conf.py.in
EXTRA_DIST = ./lofreq/__init__.py \
	./lofreq/conf.py \
	./lofreq/em.py \
	./lofreq/multiple_testing.py \
	./lofreq/qual.py \
	./lofreq/sam.py \
	./lofreq/simple_vcf.py \
	./lofreq/snp.py \
	./lofreq/utils.py \
	./lofreq_ext/cdflib90/cdflib.h \
	./lofreq_ext/cdflib90/dcdflib.c \
	./lofreq_ext/cdflib90/ipmpar.c \
	./lofreq_ext/fet.c \
	./lofreq_ext/fet.h \
	./lofreq_ext/lofreq_ext.c \
	./scripts/lofreq_alnoffset.py \
	./scripts/lofreq_bonf.py \
	./scripts/lofreq_diff.py \
	./scripts/lofreq_filter.py \
	./scripts/lofreq_pileup_summary.py \
	./scripts/lofreq_regionbed.py \
	./scripts/lofreq_snpcaller.py \
	./scripts/lofreq_uniq.py \
	./scripts/lofreq_uniq_pipeline.py \
	./scripts/lofreq_varpos_to_vcf.py \
	./scripts/lofreq_version.py \
	./setup.py \
	./setup_conf.py \
	./utils/lofreq_ign_pos.py \
	./utils/lofreq_others.py

all:
	 $(PYTHON) setup.py build
install:
	 $(PYTHON) setup.py install --prefix '$(prefix)'

