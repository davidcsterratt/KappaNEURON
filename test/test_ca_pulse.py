## HH cell with spine and excitatory autapase, along with an L-type Ca
## channel on the spine head, calcium accumulation and an
## reaction-diffusion (rxd) mechanism.
from neuron import *
from neuron import rxd

## Neuron

# Spine Head
sh = h.Section()
sh.insert("pas")                # Passive channel
sh.insert("capulse")            # Code to give Ca pulse
sh.L = 0.1
sh.diam = 4

## This setting of parameters gives a calcium influx and pump
## activation that is more-or-less scale-independent
sh.gcalbar_capulse = 0.05*sh.diam

## Reaction-diffusion mechanism
## This appears to integrate the incoming Ca
# WHERE the dynamics will take place
r = rxd.Region([sh], nrn_region='i')

# WHO are the actors
ca = rxd.Species(r, name='ca', charge=2, initial=0.01)
P  = rxd.Species(r, name='P',  charge=0, initial=0.2)
kappa = rxd.Kappa([ca, P], "caPump.ka", r)
vol = sh.L*numpy.pi*(sh.diam/2)**2
kappa.setVariable('gamma1', 1E-3*0.1*numpy.pi*0.5**2/vol)

## Current clamp stimulus
stim = h.VClamp(sh(0.5))
stim.dur[0] = 10
stim.dur[1] = 10
stim.dur[2] = 10
stim.amp[0] = -70
stim.amp[1] =  0
stim.amp[2] = -70

## Record Time from NEURON (neuron.h._ref_t)
rec_t = h.Vector()
rec_t.record(h._ref_t)
## Record Voltage from the center of the soma
rec_v = h.Vector()
rec_v.record(sh(0.5)._ref_v)
## Record Ca from spine head
rec_cai = h.Vector()
rec_cai.record(sh(0.5)._ref_cai)
## Record ica from spine head
rec_ica = h.Vector()
rec_ica.record(sh(0.5)._ref_ica)
## Record P from spine head
rec_Pi = h.Vector()
rec_Pi.record(sh(0.5)._ref_Pi)

## Run
init()
run(30)

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
ica = []
Pi = []
times.append(list(rec_t)) # alternativ to `list(rec_t)`: `numpy.array(rec_t)`
voltages.append(list(rec_v))
cai.append(list(rec_cai))
ica.append(list(rec_ica))
Pi.append(list(rec_Pi))
# check types by:
# >>> type(rec_t)
# >>> type(time[0])

fig, (ax1, ax2, ax3, ax4) = plt.subplots(nrows=4, ncols=1)
ax1.plot(times[0], voltages[0])
ax1.set_xlabel("Time [ms]")
ax1.set_ylabel("V [mV]")
ax1.axis(ymin=-80, ymax=50)

ax2.plot(times[0], ica[0])
ax2.set_xlabel("Time [ms]")
ax2.set_ylabel("ICa [mA/cm2]")
ax2.axis(ymin=-1, ymax=0.1)

ax3.plot(times[0], cai[0])
ax3.set_xlabel("Time [ms]")
ax3.set_ylabel("Ca [mM]")
ax3.axis(ymin=-1E-2, ymax=0.5E-1)

ax4.plot(times[0], Pi[0])
ax4.set_xlabel("Time [ms]")
ax4.set_ylabel("P [mM]")
ax4.axis(ymin=-1E-5, ymax=2.5E-1)

fig.show() # If the interpreter stops now: close the figure.
# For interactive plotting, see `Part 1` -> `ipython`

fig.savefig("../doc/test_ca_pulse.pdf", format='pdf')

numpy.savez("test_ca_pulse", t=times[0], cai=cai[0])
