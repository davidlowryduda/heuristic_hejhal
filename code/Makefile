

.PHONY: usage
usage:
	@echo
	@echo " MAKE USAGE:"
	@echo "   make cython --- recompile *.spyx files"
	@echo "   make python --- recompile *.py files from *.sage files"
	@echo "   make tests  --- run tests"
	@echo

%.cpython-38-x86_64-linux-gnu.so: %.spyx
	sage -c "from sage.misc.cython import cython; cython('$<', compile_message=True, use_cache=False, create_local_so_file=True)"

.PHONY: cython
cython: lpkbessel.cpython-38-x86_64-linux-gnu.so \
	      maass_levelone_computations.cpython-38-x86_64-linux-gnu.so \
				maass_levelone_computations_tests.cpython-38-x86_64-linux-gnu.so \
				maass_sqfreelevel_computations.cpython-38-x86_64-linux-gnu.so \
				maass_sqfreelevel_tests.cpython-38-x86_64-linux-gnu.so

%.py: %.sage
	sage -preparse $<
	mv $<.py $@

.PHONY: python
python: groupdata.py

.PHONY: tests
tests: python cython
	# sage levelone_test_runner.sage
	sage sqfreelevel_test_runner.sage
