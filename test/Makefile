PREFIX = /home/dcs
# PREFIX = /disk/scratch/sterratt
PYTHONPATH = $(PREFIX)/lib/python

MODFILES = *.mod *.inc

x86_64/special: $(MODFILES)
	$(PREFIX)/x86_64/bin/nrnivmodl

reference: x86_64/special
	../nrn-kappa test_ca_pulse_mod.py

test: reference
	@ ../nrn-kappa test_ca_pulse.py
