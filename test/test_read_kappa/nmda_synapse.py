## Spine head into which Ca flows via an NMDA channel and is pumped
## out. Voltage clamp ensures almost constant Ca flow.
from neuron import *
from neuron import rxd
import numpy

## Neuron
h.celsius = 34                  # Temperature
g_pas = 0.0001                  # Gives membrane time constant of 10ms

# Spine Head
sh = h.Section()
sh.insert("pas")                # Passive channel
# This gives volume of 0.1um3
sh.L = 0.2
sh.diam = 0.4
sh.g_pas = g_pas

# Synapse
synstim = h.NetStim()
synstim.start = 10
synstim.noise = 0
synstim.number = 10
synstim.interval = 25

#ampasyn     = h.AmpaSyn(sh(0.5))
#ampanetcon  = h.NetCon(synstim, ampasyn)
#ampanetcon.weight[0] = 0.20E-3     # From the ddsp work

nmdasyn     = h.NmdaSyn(sh(0.5))
h.K0_NmdaSyn =  2.57                   
h.delta_NmdaSyn = 0.96
nmdanetcon  = h.NetCon(synstim, nmdasyn)
nmdanetcon.weight[0] = 0.045E-3     # From the ddsp work

## Reaction-diffusion mechanism
## This appears to integrate the incoming Ca
# WHERE the dynamics will take place
r = rxd.Region([sh], nrn_region='i')

# WHO are the actors
# ca        = rxd.Species(r, name='ca'   , charge=2, initial=0.001)
NMDA      = rxd.Species(r, name='NMDA'  , charge=1, initial=0.01)
Glu       = rxd.Species(r, name='Glu'   , charge=1, initial=0)
#NMDAO     = rxd.Species(r, name='NMDAO' , charge=0)

kappa = rxd.Kappa([NMDA, Glu], "nmda.ka", r, time_units="ms", verbose=True)

## This setting of parameters gives a calcium influx and pump
## activation that is more-or-less scale-independent
vol = sh.L*numpy.pi*(sh.diam/2)**2

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
## Record Ca bound to CB from spine head
rec_NMDAi = h.Vector()
rec_NMDAi.record(sh(0.5)._ref_NMDAi)
## Record Free CaM from spine head
rec_Glui = h.Vector()
rec_Glui.record(sh(0.5)._ref_Glui)

## Run
init()
#kappa.run_free(200)
print("Running NEURON-kappa")
run(20)

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
Glui = []

times.append(list(rec_t)) # alternativ to `list(rec_t)`: `numpy.array(rec_t)`
voltages.append(list(rec_v))
cai.append(list(rec_cai))
ica.append(list(rec_ica))
Glui.append(list(rec_Glui))

def plot_data(tmax=None):
    if (tmax == None): 
        tmax = max(times[0])
    
    fig, (ax1, ax2, ax3, ax4) = plt.subplots(nrows=4, ncols=1)
    ax1.plot(times[0], voltages[0])
    ax1.set_xlabel("Time [ms]")
    ax1.set_ylabel("V [mV]")
    ax1.axis(ymin=-80, ymax=50)
    ax1.axis(xmin=0, xmax=tmax)

    ax2.plot(times[0], ica[0])
    ax2.set_xlabel("Time [ms]")
    ax2.set_ylabel("ICa [mA/cm2]")
    ## ax2.axis(ymin=-1, ymax=0.1)
    ax2.axis(xmin=0, xmax=tmax)

    ax3.plot(times[0], cai[0])
    ax3.plot(times[0], rec_NMDAi)
    ax3.plot(times[0], rec_Glui)
    ax3.set_xlabel("Time [ms]")
    ax3.set_ylabel("[mM]")
    plt.axes(ax3)
    plt.legend(('Ca', 'NMDA', 'Glu'))
    ## ax3.axis(ymin=-1E-2, ymax=0.5E-1)
    ax3.axis(xmin=0, xmax=tmax)

    fig.show() # If the interpreter stops now: close the figure.
    # For interactive plotting, see `Part 1` -> `ipython`

    # fig.savefig("../doc/simple-psd-pepke.pdf", format='pdf')

numpy.savez("nmda_synapse.npz", t=times[0], cai=cai[0],  ica=ica[0], voltages=voltages[0], diam=sh.diam)

