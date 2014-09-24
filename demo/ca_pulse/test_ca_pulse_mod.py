## Spine head into which Ca flows via a very simple channel and is
## pumped out.  Voltage clamp ensures almost constant Ca flow.
from test_ca_pulse_common import *
import numpy

def run(diam=0.2, 
        gbar=0.05,
        gamma2=1,
        P0=0.2,
        vclamp=False):
    # Neuron spine Head
    sh = make_spine_head(diam=diam, gbar=gbar)
    sh.insert("caPump2")            # My own calcium buffer

    ## This setting of parameters gives a calcium influx and pump
    ## activation that is more-or-less scale-independent
    vol = sh.L*numpy.pi*(sh.diam/2)**2 # Volume in um3
    ## sh.gamma1_caPump = 1E-3*(0.1*numpy.pi*((1./2)**2))/vol
    ## sh.gamma2_caPump = gamma2
    sh.k1_caPump2 = 47.3
    sh.k2_caPump2 = gamma2

    sh.P0_caPump2 = P0

    if vclamp:
        stim = insert_vclamp(sh)

    ## Record P from spine head
    rec_Pi = h.Vector()
    rec_Pi.record(sh(0.5)._ref_P_caPump2)

    run_and_save(sh, rec_Pi, 'test_ca_pulse_mod')
