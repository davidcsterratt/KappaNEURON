Demo of calcium accumulation
============================

Basic example
-------------

```
rm -Rf x86_64
nrnivmodl
python -i basic.py
```

Comparison of kappa and mod file implementation
-----------------------------------------------

To run a comparision of the kappa simulation with the equivalent mod
simulation:
```
rm -Rf x86_64
nrnivmodl
python -i test_ca_pulse_run.py
```
