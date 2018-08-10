#!/usr/bin/python2.7
import neuron
import os

## Find the version of nrnivmodl corresponding to this python module
neuron_root = os.path.realpath(os.path.join(os.path.dirname(os.path.realpath(neuron.__file__)), '../../../'))
## Use it to compile files in the test directory
os.system(os.path.join(neuron_root, 'x86_64/bin/nrnivmodl test'))
