# This is really just a wrapper to setup.py since that script takes
# care of building and installing the Python stuff
#
# Python scripts and other source has to be added manually here to
# make sure they are incorporated when calling 'make dist'
#
# How to get this list:
# source=$(python setup.py --dry-run sdist 2>/dev/null | grep 'hard linking'  | cut -f 3 -d ' ');
# source="$source setup_conf.py $(ls ./utils/*py)";
# echo $source | sed -e 's, , \\\n,g' | sed -e 's,^,\t,';
#
EXTRA_DIST = setup.py \
	lofreq/__init__.py \
	lofreq/conf.py \
	lofreq/em.py \
	lofreq/multiple_testing.py \
	lofreq/posmap.py \
	lofreq/qual.py \
	lofreq/sam.py \
	lofreq/simple_vcf.py \
	lofreq/snp.py \
	lofreq/utils.py \
	lofreq_ext/lofreq_ext.c \
	scripts/lofreq_alnoffset.py \
	scripts/lofreq_bonf.py \
	scripts/lofreq_diff.py \
	scripts/lofreq_filter.py \
	scripts/lofreq_pileup_summary.py \
	scripts/lofreq_regionbed.py \
	scripts/lofreq_snpcaller.py \
	scripts/lofreq_uniq.py \
	scripts/lofreq_uniq_pipeline.py \
	scripts/lofreq_varpos_to_vcf.py \
	scripts/lofreq_version.py \
	setup_conf.py \
	./utils/lofreq_ign_pos.py \
	./utils/lofreq_others.py

all:
	$(PYTHON) setup.py build

# make sure to actually install stuff via python's setup.py
install-exec-hook:
	$(PYTHON) setup.py install --prefix '$(prefix)'

# local clean target: call setup.py and remove ./build/
clean-local:
	$(PYTHON) setup.py clean
	rm -rf ./build
