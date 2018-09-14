from neuron import *
import KappaNEURON
from neuron import rxd
import matplotlib.pyplot as plt

## Create a compartment with a passive channel
sh = h.Section()
sh.L = 0.1
sh.diam = 0.2
sh.insert('pas')

# Code to give Ca pulse from 5ms to 10ms
sh.insert('capulse')
sh.gbar_capulse = 0.001
sh.fghk_capulse = 1         # Use GHK
sh.t0_capulse = 5
sh.t1_capulse = 10

## Define region where the dynamics will occur ('i' means intracellular)
r = rxd.Region([sh], nrn_region='i')

## Define the species, the ca ion (already built-in to NEURON), and
## the pump molecule. These names must correspond to the agent names
## in the Kappa file.
ca = rxd.Species(r, name='ca', charge=2, initial=0.0)
P  = rxd.Species(r, name='P',  charge=0, initial=0.2)

## Create the link between the Kappa model and the species just defined
kappa = KappaNEURON.Kappa(membrane_species=[ca], species=[P],
kappa_file='caPump2.ka', regions=[r])

## Transfer variable settings to the kappa model
kappa.setVariable('k1', 47.3)
kappa.setVariable('k2', 0.1)

## Volume needs to be specified explicitly
vol = sh.L*numpy.pi*(sh.diam/2)**2
kappa.setVariable('vol', vol)

## Set up recordings

## Record Time from NEURON (neuron.h._ref_t)
rec_t = h.Vector()
rec_t.record(h._ref_t)
## Record Voltage from the center of the soma
rec_v = h.Vector()
rec_v.record(sh(0.5)._ref_v)
## Record ica from spine head
rec_ica = h.Vector()
rec_ica.record(sh(0.5)._ref_ica)
## Record Ca from spine head
rec_cai = h.Vector()
rec_cai.record(sh(0.5)._ref_cai)
## Record P from spine head
rec_Pi = h.Vector()
rec_Pi.record(sh(0.5)._ref_Pi)

## Run
init()
h.finitialize(sh.e_pas)
run(30)

## Plot
fig, ax = plt.subplots(nrows=4, ncols=1)
ax[0].plot(numpy.array(rec_t), numpy.array(rec_v), 'r')
ax[1].plot(numpy.array(rec_t), numpy.array(rec_ica), 'r')
ax[2].plot(numpy.array(rec_t), numpy.array(rec_cai), 'r')
ax[3].plot(numpy.array(rec_t), numpy.array(rec_Pi), 'r')
fig.show()
