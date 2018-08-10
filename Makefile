# .PHONY test: all
#	make -C test/test_ca_pulse test

PREFIX=/usr/local

build_py:
	python2.7 setup.py build

.PHONY test: build_py
	./KappaNEURON/tests/nrnivmodl.py
	python2.7 -m unittest KappaNEURON.tests.TestCaAccumulation.test_injectCalcium && \
	python2.7 -m unittest KappaNEURON.tests.TestCaAccumulation.test_injectCalciumGHK && \
	python2.7 -m unittest KappaNEURON.tests.TestCaAccumulation.test_injectCalciumPump && \
	python2.7 -m unittest KappaNEURON.tests.TestCaAccumulation.test_injectCalciumPumpGHK && \
	python2.7 -m unittest KappaNEURON.tests.TestCaAccumulation.test_injectCalciumPump2 && \
	python2.7 -m unittest KappaNEURON.tests.TestCaAccumulation.test_injectCalciumPump2k2 && \
	python2.7 -m unittest KappaNEURON.tests.TestCaAccumulation.test_twoMembraneSpecies && \
	python2.7 -m unittest KappaNEURON.tests.TestCaAccumulation.test_twoMembraneSpeciesOneUncharged && \
	echo "All tests passed"

install:
	python2.7 setup.py install --prefix=$(PREFIX)

