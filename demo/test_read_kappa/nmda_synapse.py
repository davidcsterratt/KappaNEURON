## Spine head into which Ca flows via an NMDA channel and is pumped
## out. Voltage clamp ensures almost constant Ca flow.
from neuron import *
from neuron import rxd
import numpy

## Neuron
h.celsius = 34                  # Temperature
g_pas = 0.0001                  # Gives membrane time constant of 10ms

# Dendrite
dend = h.Section()
dend.insert("pas")                # Passive channel
dend.L = 100
dend.diam = 1
dend.g_pas = g_pas
dend.insert("hh")

# Spine neck needed probably!
sn = h.Section()
sn.insert("pas")                # Passive channel
sn.L = 1
sn.diam = 0.1
sn.g_pas = g_pas
sn.connect(dend, 0.5, 0)

# Spine Head
sh = h.Section()
sh.insert("pas")                # Passive channel
# This gives volume of 0.1um3
sh.L = 0.2
sh.diam = 0.8
sh.g_pas = g_pas
sh.connect(sn, 1, 0)

# Synapse
synstim = h.NetStim()
synstim.start = 10
synstim.noise = 0
synstim.number = 10
synstim.interval = 25

ampasyn     = h.AmpaSyn(sh(0.5))
ampanetcon  = h.NetCon(synstim, ampasyn)
ampanetcon.weight[0] = 0.20E-3     # From the ddsp work

nmdasyn     = h.NmdaSyn(sh(0.5))
h.K0_NmdaSyn =  2.57                   
h.delta_NmdaSyn = 0.96
nmdanetcon  = h.NetCon(synstim, nmdasyn)
nmdanetcon.weight[0] = 0.045E-3     # From the ddsp work

## This setting of parameters gives a calcium influx and pump
## activation that is more-or-less scale-independent
vol = sh.L*numpy.pi*(sh.diam/2)**2
N_A = 6.02205E23 # Avogadro's constant
# Concentration of one agent in the volume in mM 
agconc = 1E18/(N_A * vol)
## Reaction-diffusion mechanism
## This appears to integrate the incoming Ca
# WHERE the dynamics will take place
r = rxd.Region([sh], nrn_region='i')

# WHO are the actors
ca        = rxd.Species(r, name='ca'       , charge=2, initial=0.001)
NMDA      = rxd.Species(r, name='NMDA'     , charge=0, initial=3*agconc)
NMDAC0    = rxd.Species(r, name='NMDAC0' , charge=0)
NMDAC1    = rxd.Species(r, name='NMDAC1' , charge=0)
NMDAC2    = rxd.Species(r, name='NMDAC2' , charge=0)
NMDAC3    = rxd.Species(r, name='NMDAC3' , charge=0)
Glu       = rxd.Species(r, name='Glu'   , charge=1, initial=0)
#NMDAO     = rxd.Species(r, name='NMDAO' , charge=0)

kappa = rxd.Kappa([NMDA, Glu], "nmda2.ka", r, time_units="ms") #, verbose=True)
# rxd.rxd.verbose=False


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
rec_NMDAC0i = h.Vector()
rec_NMDAC0i.record(sh(0.5)._ref_NMDAC0i)
rec_NMDAC1i = h.Vector()
rec_NMDAC1i.record(sh(0.5)._ref_NMDAC1i)
rec_NMDAC2i = h.Vector()
rec_NMDAC2i.record(sh(0.5)._ref_NMDAC2i)
rec_NMDAC3i = h.Vector()
rec_NMDAC3i.record(sh(0.5)._ref_NMDAC3i)


## Record Free CaM from spine head
rec_Glui = h.Vector()
rec_Glui.record(sh(0.5)._ref_Glui)
## Record Free CaM from spine head
rec_iGlu = h.Vector()
rec_iGlu.record(sh(0.5)._ref_iGlu)
## Record Free CaM from spine head
rec_iNMDA = h.Vector()
rec_iNMDA.record(sh(0.5)._ref_iNMDA)

rec_g = h.Vector()
rec_g.record(nmdasyn._ref_g)


## Run
init()
#kappa.run_free(200)
print("Running NEURON-kappa")
run(600)

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
iGlu = []
iNMDA = []

times.append(list(rec_t)) # alternativ to `list(rec_t)`: `numpy.array(rec_t)`
voltages.append(list(rec_v))
cai.append(list(rec_cai))
ica.append(list(rec_ica))



def plot_data(tmax=None):
    if (tmax == None): 
        tmax = max(times[0])
    
    fig, (ax1, ax2, ax3, ax4) = plt.subplots(nrows=4, ncols=1)
    ax1.plot(times[0], voltages[0])
    ax1.set_xlabel("Time [ms]")
    ax1.set_ylabel("V [mV]")
    ax1.axis(ymin=-66, ymax=-62)
    ax1.axis(xmin=0, xmax=tmax)

    ax2.plot(times[0], ica[0])
    ax2.plot(times[0], rec_iNMDA)
    ax2.plot(times[0], rec_iGlu)
    plt.axes(ax2)
    plt.legend(('Ca', 'NMDA', 'Glu'))
    ax2.set_xlabel("Time [ms]")
    ax2.set_ylabel("ICa [mA/cm2]")
    ## ax2.axis(ymin=-1, ymax=0.1)
    ax2.axis(xmin=0, xmax=tmax)

    ax3.plot(times[0], rec_cai)
    ax3.plot(times[0], rec_NMDAi)
    ax3.plot(times[0], rec_Glui)
    ax3.set_xlabel("Time [ms]")
    ax3.set_ylabel("[mM]")
    plt.axes(ax3)
    plt.legend(('Ca', 'NMDA', 'Glu'))
    ax3.axis(ymin=0, ymax=0.1)
    ax3.axis(xmin=0, xmax=tmax)

    ax4.plot(times[0], rec_g)
    ax4.set_xlabel("Time [ms]")
    ax4.set_ylabel("g [uS]")
    ## ax3.axis(ymin=-1E-2, ymax=0.5E-1)
    ax4.axis(xmin=0, xmax=tmax)


    fig.show() # If the interpreter stops now: close the figure.
    # For interactive plotting, see `Part 1` -> `ipython`

    # fig.savefig("../doc/simple-psd-pepke.pdf", format='pdf')

numpy.savez("nmda_synapse.npz", t=times[0], cai=cai[0],  ica=ica[0], voltages=voltages[0], diam=sh.diam)

plot_data()
