## HH cell with spine and excitatory autapase, along with an L-type Ca
## channel on the spine head, calcium accumulation and an
## reaction-diffusion (rxd) mechanism.
from neuron import *
from neuron import rxd
## Load mechanisms
## load_mechanisms('.')



## Neuron
soma = h.Section()
soma.insert("hh")

## Spine
# Neck
sn = h.Section()
sn.insert("pas")
sn.L = 2
sn.diam = 0.1

# Head
sh = h.Section()
sh.insert("pas")                # Passive channel
## sh.insert("cacum")              # Ca accumulation
sh.insert("cal")                # L-type Ca channel
sh.gcalbar_cal = 1
sh.L = 0.1
sh.diam = 1

## Connections
sn.connect(soma, 0.5, 0) # Connect soma(0.5) with sn(0)
sh.connect(sn, 1, 0) # Connect sh(1) with sh(0)

## Autapse
aut = h.ExpSyn(sn(1))
nc = h.NetCon(soma(0.5)._ref_v, aut, -20, 3, 10)

## Reaction-diffusion mechanism
## This appears to integrate the incoming Ca
# WHERE the dynamics will take place
r = rxd.Region([sh], nrn_region='i')

# WHO are the actors
ca = rxd.Species(r, name='ca', charge=2, initial=0.01)

## rate = rxd.Rate2(ca, -0.1*ca)

## HOW do they react
## reaction = rxd.Rection(ca > , 1)

## Current clamp stimulus
stim = h.IClamp(soma(0.5))
stim.delay = 1
stim.amp = 10
stim.dur = 3

## Record Time from NEURON (neuron.h._ref_t)
rec_t = h.Vector()
rec_t.record(h._ref_t)
## Record Voltage from the center of the soma
rec_v = h.Vector()
rec_v.record(soma(0.5)._ref_v)
## Record Ca from spine head
rec_cai = h.Vector()
rec_cai.record(sh(0.5)._ref_cai)

## Run
init()
run(10)

## Plot
import matplotlib.pyplot as plt

## Get values from NEURON-vector format into Python format
times = [] # Use list to add another trace later.
voltages = []
times.append(list(rec_t)) # alternative to `list(rec_t)`:
                          # `numpy.array(rec_t)`
voltages.append(list(rec_v))

# Plot the recordings with matplotlib
# ===================================

import matplotlib.pyplot as plt

# get values from NEURON-vector format into Python format
times = [] # Use list to add another trace later.
voltages = []
cai = []
times.append(list(rec_t)) # alternativ to `list(rec_t)`: `numpy.array(rec_t)`
voltages.append(list(rec_v))
cai.append(list(rec_cai))
# check types by:
# >>> type(rec_t)
# >>> type(time[0])

fig, (ax1, ax2) = plt.subplots(nrows=2, ncols=1)
ax1.plot(times[0], voltages[0])
ax1.set_xlabel("Time [ms]")
ax1.set_ylabel("Voltage [mV]")
ax1.axis(ymin=-80, ymax=50)

ax2.plot(times[0], cai[0])
ax2.set_xlabel("Time [ms]")
ax2.set_ylabel("Ca [mM]")
ax2.axis(ymin=0, ymax=0.5)


fig.show() # If the interpreter stops now: close the figure.
# For interactive plotting, see `Part 1` -> `ipython`




