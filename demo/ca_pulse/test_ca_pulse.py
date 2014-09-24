## Spine head into which Ca flows via a very simple channel and is
## pumped out.  Voltage clamp ensures almost constant Ca flow.
from test_ca_pulse_common import *
from neuron import rxd
import KappaNEURON

def run(diam=0.2, 
        gbar=0.05,
        gamma2=1,
        P0=0.2,
        vclamp=False):
    # Make spine head in NEURON
    sh = make_spine_head(diam=diam, gbar=gbar)

    ## Create region where the dynamics will take place
    r = rxd.Region([sh], nrn_region='i')

    ## Create species
    # Calcium
    ca = rxd.Species(r, name='ca', charge=2, initial=0.0)

    # Pump molecule
    P  = rxd.Species(r, name='P',  charge=0, initial=P0)

    ## Create Kappa simulation defined in caPump.ka in the context of
    ## the spine. Since Calcium crosses the membrane it is given in
    ## the membrane_species argument, whereas the pump molecule is
    ## defined in the species argument as it is purely internal. 
    KappaNEURON.mode='continuous_influx'
    kappa = KappaNEURON.Kappa(membrane_species=[ca], species=[P], kappa_file='caPump2.ka', regions=r)
    ## Set variables in the Kappa simulation
    vol = sh.L*numpy.pi*(sh.diam/2)**2
    kappa.setVariable('k1', 47.3)
    kappa.setVariable('k2', gamma2)
    kappa.setVariable('vol', vol)

    ## Insert Vclamp if desired
    if vclamp:
        stim = insert_vclamp(sh)

    ## Record P from spine head
    rec_Pi = h.Vector()
    rec_Pi.record(sh(0.5)._ref_Pi)
    
    run_and_save(sh, rec_Pi, 'test_ca_pulse')
