x86_64/special: capulse.mod caPump.mod
	cd test; /disk/scratch/sterratt/x86_64/bin/nrnivmodl

reference:
	python2.7 test_ca_pulse_mod.py

test: reference
	@ CLASSPATH=..:/home/sterratt/projects/mlm/SpatialKappa/SpatialKappa-v2.1.1.jar:'../../SpatialKappa/*':../../py4j/py4j0.7.jar python2.7 -i test_ca_pulse.py
