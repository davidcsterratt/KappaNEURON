## Spine head into which Ca flows via an NMDA channel and is pumped
## out. Voltage clamp ensures almost constant Ca flow.
from neuron import *
from neuron import rxd
import KappaNEURON
import numpy
import random
import os

## Create Neuron
h.celsius = 34                  # Temperature
g_pas = 0.0001               # Gives membrane time constant of 10ms

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
#h.dt = 0.001
sh.g_pas = g_pas
sh.connect(sn, 1, 0)

## Create stimuluse to trigger synapse containing Kappa mechanism
## paired with action potential
synstim = h.NetStim()
synstim.start = 100
synstim.noise = 0
synstim.number = 10
synstim.interval = 50

## Create the AMPA and NMDA synapses
ampasyn     = h.AmpaSyn(sh(0.5))
ampanetcon  = h.NetCon(synstim, ampasyn)
ampanetcon.weight[0] = 0.20E-3

nmdasyn     = h.NmdaSynUrak(sh(0.5))
h.K0_NmdaSynUrak =  2.57                   
h.delta_NmdaSynUrak = 0.96
nmdanetcon  = h.NetCon(synstim, nmdasyn)
nmdanetcon.weight[0] = 0.045E-3

## Used to initiate action potential 1ms after the synaptic stimulus
apinit = h.ExpSyn(dend(0))
apinit.tau = 0.3/3
apinit.e = 0
apinitcon = h.NetCon(synstim, apinit)
apinitcon.weight[0] = 3*0.4/70 # 0.4nA/70mV = 0.4/70 uS
apinitcon.delay = 1

## Create synaptic bombardment
class MyAmpaSyn():
    def __init__(self):
        self.stim = h.NetStim()
        self.stim.start = 20
        self.stim.noise = 1
        self.stim.number = 400
        self.stim.interval = 20
        
        self.syn = h.AmpaSyn(dend(random.random()))
        self.netcon  = h.NetCon(self.stim, self.syn)
        self.netcon.weight[0] = 0.20E-3     # From the ddsp work

synlist = []
for i in range(1, 50):
    synlist.append(MyAmpaSyn())

## This setting of parameters gives a calcium influx and pump
## activation that is more-or-less scale-independent
vol = sh.L*numpy.pi*(sh.diam/2)**2
N_A = 6.02205E23 # Avogadro's constant
# Concentration of one agent in the volume in mM 
agconc = 1E18/(N_A * vol)

## Create region where the dynamics will take place
r = rxd.Region([sh], nrn_region='i')

## Create species
ca         = rxd.Species(r, name='ca'        , charge=2, initial=0.001)
Glu        = KappaNEURON.UnchargedSpecies(r, name='Glu', initial=0)
NMDA       = rxd.Species(r, name='NMDA'      , charge=0, initial=19*agconc)
CB         = rxd.Species(r, name='CB'        , charge=0, initial=0.100) # Faas &al
cam        = rxd.Species(r, name='CaM'       , charge=0, initial=0.030) # Fass &al, Pepke &al
CaMKII     = rxd.Species(r, name='CaMKII'    , charge=0, initial=0.080) # Pepke &al
CaCB       = rxd.Species(r, name='CaCB'      , charge=0) 
CaCaMC     = rxd.Species(r, name='CaCaMC'    , charge=0)
CaCaMN     = rxd.Species(r, name='CaCaMN'    , charge=0)
KCaCaM2C   = rxd.Species(r, name='KCaCaM2C'  , charge=0)
CaMKIIp    = rxd.Species(r, name='CaMKIIp'   , charge=0)
stargazinp = rxd.Species(r, name='stargazinp', charge=0)

## Create Kappa simulation defined in caPump.ka in the context of
## the spine. Since Calcium crosses the membrane it is given in
## the membrane_species argument, whereas the pump molecule is
## defined in the species argument as it is purely internal. 
kappa = KappaNEURON.Kappa(membrane_species=[ca, Glu], species=[NMDA, CB, cam, CaMKII, CaCB, CaCaMC, CaCaMN, KCaCaM2C, CaMKIIp, stargazinp], kappa_file='simple-psd-pepke-kappa-nmda.ka', regions=r)
rxd.rxd.verbose=False

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
## Record iGlu from spine head
rec_iGlu = h.Vector()
rec_iGlu.record(sh(0.5)._ref_iGlu)
## Record Ca bound to CB from spine head
rec_CaCBi = h.Vector()
rec_CaCBi.record(sh(0.5)._ref_CaCBi)
## Record Free CaM from spine head
rec_cami = h.Vector()
rec_cami.record(sh(0.5)._ref_CaMi)

rec_CaMKIIi = h.Vector()
rec_CaMKIIi.record(sh(0.5)._ref_CaMKIIi)
rec_CaCaMNi = h.Vector()
rec_CaCaMNi.record(sh(0.5)._ref_CaCaMNi)
rec_CaCaMCi = h.Vector()
rec_CaCaMCi.record(sh(0.5)._ref_CaCaMCi)
rec_KCaCaM2Ci = h.Vector()
rec_KCaCaM2Ci.record(sh(0.5)._ref_KCaCaM2Ci)
rec_CaMKIIpi = h.Vector()
rec_CaMKIIpi.record(sh(0.5)._ref_CaMKIIpi)
rec_stargazinpi = h.Vector()
rec_stargazinpi.record(sh(0.5)._ref_stargazinpi)


## Run
init()
print("Running kappa-only to initialise")
## FIXME: put in some read-out to check when system has equilibriated.
kappa.run_free(100)
print("Running NEURON-kappa")
run(6000)
if (0):
    for i in range(1,60):
        print("Running kappa-only")
        kappa.run_free(990)
        print("Running NEURON-kappa and trying to trick NEURON")
        h.t = h.t + 990
        run(h.t + 10)

# Plot the recordings with matplotlib
# ===================================

import matplotlib.pyplot as plt

# Get values from NEURON-vector format into Python format
times = [] 
voltages = []
cai = []
ica = []
cami = []

times.append(list(rec_t)) # alternative to `list(rec_t)`: `numpy.array(rec_t)`
voltages.append(list(rec_v))
cai.append(list(rec_cai))
ica.append(list(rec_ica))
cami.append(list(rec_cami))

def plot_data(tmax=None):
    if (tmax == None): 
        tmax = max(times[0])
    
    fig, (ax1, ax2, ax3, ax4, ax5) = plt.subplots(nrows=5, ncols=1)
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

    ax3.plot(times[0], rec_cai)
    ax3.plot(times[0], rec_CaCBi)
    ax3.plot(times[0], rec_CaCaMNi)
    ax3.plot(times[0], rec_CaCaMCi)
    ax3.set_xlabel("Time [ms]")
    ax3.set_ylabel("[mM]")
    plt.axes(ax3)
    plt.legend(('Ca', 'CaCB', 'CaCaMN', 'CaCaMC'))
    ax3.axis(ymin=-1E-3, ymax=0.5E-1)
    ## ax3.axis(ymin=-1E-2, ymax=0.5E-1)
    ax3.axis(xmin=0, xmax=tmax)

    
    ax4.plot(times[0], rec_KCaCaM2Ci)
    ax4.plot(times[0], rec_CaMKIIpi)
    ax4.plot(times[0], rec_stargazinpi)
    ax4.set_xlabel("Time [ms]")
    ax4.set_ylabel("[mM]")
    plt.axes(ax4)
    plt.legend(('KCaCaM2C','CaMKIIp', 'stargazinp'))
    ax4.axis(xmin=0, xmax=tmax)
    ax4.axis(ymin=0, ymax=0.001)

    ax5.plot(times[0], numpy.array(rec_iGlu))

    fig.show() # If the interpreter stops now: close the figure.
    # For interactive plotting, see `Part 1` -> `ipython`
    try:
        os.makedirs("figs")
    except Exception as e:
        pass
        
    fig.savefig("figs/simple-psd-pepke-kappa-nmda.pdf", format='pdf')

# numpy.savez("simple-psd-pepke-kappa-nmda", t=times[0], cai=rec_cai, 
#             cami=cami[0], ica=ica[0], voltages=voltages[0], 
#             diam=sh.diam, Glu=numpy.array(Glu), 
#             NMDA=numpy.array(NMDA), CB=numpy.array(CB),
#             cam=numpy.array(cam), CaMKII=numpy.array(CaMKII), 
#             CaCB=numpy.array(CaCB), CaCaMC=numpy.array(CaCaMC), 
#             CaCaMN=numpy.array(CaCaMN), KCaCaM2C=numpy.array(KCaCaM2C),
#             CaMKIIp=numpy.array(CaMKIIp), stargazinp=numpy.array(stargazinp))

numpy.savez("simple-psd-pepke-kappa-nmda-comp", times=times, rec_cai=rec_cai, 
            cami=cami, ica=ica, voltages=voltages, 
            diam=sh.diam, rec_CaCBi=rec_CaCBi,
            rec_CaMKIIi=rec_CaMKIIi, rec_CaCaMNi=rec_CaCaMNi,
            rec_CaCaMCi=rec_CaCaMCi, rec_KCaCaM2Ci=rec_KCaCaM2Ci,
            rec_CaMKIIpi=rec_CaMKIIpi,  rec_stargazinpi=rec_stargazinpi)

plot_data()

# Total Ca influx is ica * dt * area /2/F *NA

# Area = pi*diam*L = pi*0.8*0.2 um2 = pi*0.8*0.2*10^-8cm2

# ica is measured in mA/cm2; dt in ms, so the total charge in Coulombs is
# Q = 10^-3*ica * 10^-3*dt * pi*diam*L*10^-8
