PREFIX = /home/dcs
# PREFIX = /disk/scratch/sterratt

MODFILES = *.mod

x86_64/special: $(MODFILES)
	$(PREFIX)/x86_64/bin/nrnivmodl

test: x86_64/special
	PYTHONPATH=$(PREFIX)/lib/python python2.7 -i test_ca_pulse_common.py

