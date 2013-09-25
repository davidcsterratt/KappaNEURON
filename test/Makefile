PREFIX = /home/dcs
# PREFIX = /disk/scratch/sterratt
PYTHONPATH = $(PREFIX)/lib/python

x86_64/special: capulse.mod caPump.mod
	$(PREFIX)/x86_64/bin/nrnivmodl

reference: x86_64/special
	python2.7 test_ca_pulse_mod.py

test: reference
	@ PYTHONPATH=$(PREFIX)/lib/python CLASSPATH=..:$(HOME)/projects/mlm/SpatialKappa/SpatialKappa-v2.1.1.jar:'../../SpatialKappa/*':../../py4j/py4j0.7.jar python2.7 -i test_ca_pulse.py
