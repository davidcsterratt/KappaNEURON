# .PHONY test: all
#	make -C test/test_ca_pulse test

PREFIX=/usr/local

build_py:
	python2.7 setup.py build

.PHONY test: build_py
	python2.7 -m unittest test

install:
	python2.7 setup.py install --prefix=$(PREFIX)



