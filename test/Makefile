PREFIX = /home/dcs
# PREFIX = /disk/scratch/sterratt
PYTHONPATH = $(PREFIX)/lib/python

MODFILES = *.mod *.inc

x86_64/special: $(MODFILES)
	$(PREFIX)/x86_64/bin/nrnivmodl

refsim1: x86_64/special
	../nrn-kappa test_ca_pulse_mod.py
	@ ../nrn-kappa test_ca_pulse.py

test: refsim1


