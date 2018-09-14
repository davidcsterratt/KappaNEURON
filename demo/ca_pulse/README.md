Demo of calcium accumulation
============================

Basic example
-------------

The file `caPump2.ka` contains the kappa code for a simple pump
mechanism shown in Fig. 3 of [Sterratt & al (2015)][1].
It is incorporated in a model specified in `basic.py`, which is an
expanded version of the code in Fig. 4 of
[Sterratt & al, 2015][1]. To run the code, type the
following at the command line:

```
rm -Rf x86_64
nrnivmodl
python -i
basic.py
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
[1]: http://arxiv.org/abs/1411.4980 "Arxiv version of
Sterratt & al (2015)"
