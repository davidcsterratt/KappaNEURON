## Spine head into which Ca flows via an NMDA channel and is pumped
## out. Voltage clamp ensures almost constant Ca flow.
from neuron import *
from neuron import rxd
import numpy
import KappaNEURON

## Time for equilibriation in seconds
t_equil = 1

## Neuron

g_pas = 0.0001               # Gives membrane time constant of 10ms

# Dendrite
dend = h.Section()
dend.insert("pas")                # Passive channel
dend.L = 100
dend.diam = 1
dend.g_pas = g_pas
dend.insert("hh")

# Spine Head
sh = h.Section()
sh.insert("pas")                # Passive channel
sh.L = 0.1
sh.diam = 0.2
#h.dt = 0.001
sh.g_pas = g_pas
sh.connect(dend, 0.5, 0)

# Synapse
syn     = h.NmdaSyn(sh(0.5))
synstim = h.NetStim()
h.celsius = 34
synstim.start = 1001
#synstim.start = t_equil*1000 + 10
synstim.noise = 0
synstim.number = 10
synstim.interval = 25
netcon  = h.NetCon(synstim, syn)
netcon.weight[0] = 0.045E-3     # From the ddsp work
## netcon.weight[0] = 0     # From the ddsp work
h.K0_NmdaSyn =  2.57                   
h.delta_NmdaSyn = 0.96
## netcon.weight[0] = 0

## Used to initiate action potential
apinit = h.ExpSyn(dend(0))
apinit.tau = 0.3/3
apinit.e = 0
apinitcon = h.NetCon(synstim, apinit)
apinitcon.weight[0] = 3*0.4/70 # 0.4nA/70mV = 0.4/70 uS
apinitcon.delay = 16
# iclamp = h.IClamp(0.5, sec=dend)
#iclamp.delay = 17 
# iclamp.amp = 0.4 # nA
# iclamp.dur = 0.3
# Hence 0.3*0.4 = 0.12pC of charge

## Reaction-diffusion mechanism
## This appears to integrate the incoming Ca
# WHERE the dynamics will take place
r = rxd.Region([sh], nrn_region='i')

# WHO are the actors
ca         = rxd.Species(r, name='ca'         , charge=2, initial=0.001)
cam        = rxd.Species(r, name='CaM'        , charge=0, initial=0.5) # Check this
PSD95NR2   = rxd.Species(r, name='PSD95NR2'   , charge=0)
CaMKII_CaM = rxd.Species(r, name='CaMKII_CaM' , charge=0)
stargazin  = rxd.Species(r, name='stargazin'  , charge=0)
SAP97GluR1 = rxd.Species(r, name='SAP97GluR1' , charge=0)
SAP97NR2   = rxd.Species(r, name='SAP97NR2'   , charge=0)

kappa = KappaNEURON.Kappa([ca, cam, PSD95NR2, CaMKII_CaM, stargazin, SAP97GluR1, SAP97NR2], "simplx-demo.ka", r, time_units="s")


## This setting of parameters gives a calcium influx and pump
## activation that is more-or-less scale-independent
vol = sh.L*numpy.pi*(sh.diam/2)**2

## Reaction-diffusion mechanism
## This appears to integrate the incoming Ca
# WHERE the dynamics will take place
#r = rxd.Region([sh], nrn_region='i')

# WHO are the actors
#ca = rxd.Species(r, name='ca', charge=2, initial=0.01)

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
rec_cami = h.Vector()
rec_cami.record(sh(0.5)._ref_CaMi)

rec_stargazini = h.Vector()
rec_stargazini.record(sh(0.5)._ref_stargazini)
rec_SAP97GluR1i = h.Vector()
rec_SAP97GluR1i.record(sh(0.5)._ref_SAP97GluR1i)
rec_PSD95NR2i = h.Vector()
rec_PSD95NR2i.record(sh(0.5)._ref_PSD95NR2i)
rec_CaMKII_CaMi = h.Vector()
rec_CaMKII_CaMi.record(sh(0.5)._ref_CaMKII_CaMi)

## Run
init()
print("Running kappa-only to initialise")
## FIXME: put in some read-out to check when system has equilibriated.
kappa.run_free(120*1000)
# for i in range(1,t_equil):
#     run(h.t + 10)
#     kappa.run_free(990)
#     h.t = h.t + 990
print("Running NEURON-kappa")
run(3000)
for i in range(1,60):
    print("Running kappa-only")
    kappa.run_free(990)
    print("Running NEURON-kappa and trying to trick NEURON")
    h.t = h.t + 990
    run(h.t + 10)

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
cami = []

times.append(list(rec_t)) # alternativ to `list(rec_t)`: `numpy.array(rec_t)`
voltages.append(list(rec_v))
cai.append(list(rec_cai))
ica.append(list(rec_ica))
cami.append(list(rec_cami))

# check types by:
# >>> type(rec_t)
# >>> type(time[0])

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
    ax3.plot(times[0], cami[0])
    ax3.plot(times[0], rec_CaMKII_CaMi)
    ax3.set_xlabel("Time [ms]")
    ax3.set_ylabel("[mM]")
    plt.axes(ax3)
    plt.legend(('Ca', 'CaM', 'CaMKII_CaM'))
    ## ax3.axis(ymin=-1E-2, ymax=0.5E-1)
    ax3.axis(xmin=0, xmax=tmax)

    ax4.plot(times[0], rec_SAP97GluR1i)
    ax4.plot(times[0], rec_stargazini)
    ax4.set_xlabel("Time [ms]")
    ax4.set_ylabel("[mM]")
    plt.axes(ax4)
    plt.legend(('SAP97GluR2', 'stargazin'))
    ax4.axis(ymin=-1E-5, ymax=1E-2)
    ax4.axis(xmin=0, xmax=tmax)

    fig.show() # If the interpreter stops now: close the figure.
    # For interactive plotting, see `Part 1` -> `ipython`

    fig.savefig("../doc/syn_simplx.pdf", format='pdf')

numpy.savez("syn_simplx", t=times[0], cai=cai[0], cami=cami[0], ica=ica[0], voltages=voltages[0], diam=sh.diam)

