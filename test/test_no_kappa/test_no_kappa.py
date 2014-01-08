## Minimal example to test effect of KappaNEURON with no Kappa model added
import KappaNEURON
from neuron import *
import matplotlib.pyplot as plt

## Define neuron
sh = h.Section()
sh.insert("pas")                # Passive channel

## Recordings
rec_t = h.Vector()
rec_t.record(h._ref_t)
rec_v = h.Vector()
rec_v.record(sh(0.5)._ref_v)

## Set V=0 and run
sh.v = 0
init()
run(10)

## Show plot; should be decaying exponential
plt.plot(list(rec_t), list(rec_v))
plt.show()
