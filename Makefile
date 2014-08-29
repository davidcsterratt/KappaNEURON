# .PHONY test: all
#	make -C test/test_ca_pulse test

PREFIX=/usr/local

build_py:
	python2.7 setup.py build

.PHONY test: build_py
	./test/nrnivmodl.py
	python2.7 -m unittest test.TestCaAccumulation.test_injectCalcium && \
	python2.7 -m unittest test.TestCaAccumulation.test_injectCalciumGHK && \
	python2.7 -m unittest test.TestCaAccumulation.test_injectCalciumPump && \
	python2.7 -m unittest test.TestCaAccumulation.test_injectCalciumPumpGHK && \
	echo "All tests passed"

install:
	python2.7 setup.py install --prefix=$(PREFIX)

