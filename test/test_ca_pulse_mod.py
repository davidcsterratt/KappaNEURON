## Spine head into which Ca flows via a very simple channel and is
## pumped out.  Voltage clamp ensures almost constant Ca flow.
from test_ca_pulse_common import *
import numpy

def run(diam=0.2):
    # Neuron spine Head
    sh = make_spine_head(diam=diam)
    sh.insert("caPump")            # My own calcium buffer

    ## This setting of parameters gives a calcium influx and pump
    ## activation that is more-or-less scale-independent
    vol = sh.L*numpy.pi*(sh.diam/2)**2
    sh.gamma1_caPump = 1E-3*(0.1*numpy.pi*((1./2)**2))/vol
    sh.gamma2_caPump = 1

    stim = insert_vclamp(sh)

    ## Record P from spine head
    rec_Pi = h.Vector()
    rec_Pi.record(sh(0.5)._ref_P_caPump)

    run_and_save(sh, rec_Pi, 'test_ca_pulse_mod')
