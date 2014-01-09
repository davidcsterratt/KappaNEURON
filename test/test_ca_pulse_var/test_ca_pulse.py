## Spine head into which Ca flows via a very simple channel and is
## pumped out.  Voltage clamp ensures almost constant Ca flow.
from test_ca_pulse_common import *
from neuron import rxd
import KappaNEURON

def run(diam=0.2, 
        gcalbar=0.05,
        gamma2=1,
        P0=0.2):
    # Neuron spine head
    sh = make_spine_head(diam=diam, gcalbar=gcalbar)

    ## Reaction-diffusion mechanism
    ## This appears to integrate the incoming Ca
    # WHERE the dynamics will take place
    r = rxd.Region([sh], nrn_region='i')

    # WHO are the actors
    ca = rxd.Species(r, name='ca', charge=2, initial=0.0)
    P  = rxd.Species(r, name='P',  charge=0, initial=P0)
    kappa = KappaNEURON.Kappa([ca, P], "caPump.ka", r, verbose=True)
    vol = sh.L*numpy.pi*(sh.diam/2)**2
    ## kappa.setVariable('gamma1', 1E-3*0.1*numpy.pi*0.5**2/vol)
    ## kappa.setVariable('gamma2', gamma2)
    kappa.setVariable('k1', 47.3)
    kappa.setVariable('k2', gamma2)
    kappa.setVariable('vol', vol)

    stim = insert_vclamp(sh)

    ## Record P from spine head
    rec_Pi = h.Vector()
    rec_Pi.record(sh(0.5)._ref_Pi)

    run_and_save(sh, rec_Pi, 'test_ca_pulse')
