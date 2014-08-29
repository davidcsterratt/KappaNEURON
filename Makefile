# .PHONY test: all
#	make -C test/test_ca_pulse test

PREFIX=/usr/local

build_py:
	python2.7 setup.py build

.PHONY test: build_py
	./test/nrnivmodl.py
	python2.7 -m unittest test.TestCaAccumulation.test_injectCalcium
	python2.7 -m unittest test.TestCaAccumulation.test_injectCalciumGHK


install:
	python2.7 setup.py install --prefix=$(PREFIX)



