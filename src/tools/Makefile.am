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
	lofreq_star/__init__.py \
	lofreq_star/utils.py \
	scripts/lofreq2_indel_ovlp.py \
	scripts/lofreq2_local.py \
	scripts/lofreq2_vcfplot.py \
	setup_conf.py.in

all:
	$(PYTHON) setup.py build

# make sure to actually install stuff via python's setup.py
install-exec-hook:
	$(PYTHON) setup.py install --prefix '$(prefix)'
# FIXME tell users where stuff was installed and how to set PYTHONPATH

# local clean target: call setup.py and remove ./build/
clean-local:
	$(PYTHON) setup.py clean
	rm -rf ./build
