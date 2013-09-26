## Spine head into which Ca flows via a very simple channel and is
## pumped out.  Voltage clamp ensures almost constant Ca flow.
from neuron import *
from subprocess import call
import re
import matplotlib.pyplot as plt


def make_spine_head(L=0.1, diam=0.2, gcalbar=0.05):
    # Spine Head
    sh = h.Section()
    sh.insert("pas")                # Passive channel
    sh.insert("capulse")            # Code to give Ca pulse
    sh.L = L
    sh.diam = diam

    ## This setting of parameters gives a calcium influx and pump
    ## activation that is more-or-less scale-independent
    sh.gcalbar_capulse = gcalbar*sh.diam

    return(sh)

def insert_vclamp(sh):
    ## Voltage clamp stimulus
    stim = h.VClamp(sh(0.5))
    stim.dur[0] = 10
    stim.dur[1] = 10
    stim.dur[2] = 10
    stim.amp[0] = -70
    stim.amp[1] =  0
    stim.amp[2] = -70
    return(stim)

def run_and_save(sh, rec_Pi, dataname='test_ca_pulse_mod'):
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

    ## Run
    init()
    run(30)

    # get values from NEURON-vector format into Python format
    times = [] # Use list to add another trace later.
    voltages = []
    cai = []
    ica = []
    Pi = []
    times.append(list(rec_t)) # alternative to `list(rec_t)`: `numpy.array(rec_t)`
    voltages.append(list(rec_v))
    cai.append(list(rec_cai))
    ica.append(list(rec_ica))
    Pi.append(list(rec_Pi))

    numpy.savez(dataname, t=times[0], cai=cai[0], Pi=Pi[0], ica=ica[0], voltages=voltages[0], diam=sh.diam)


def compare_traces():
    import test_ca_pulse_mod
    test_ca_pulse_mod.run()
    import test_ca_pulse
    test_ca_pulse.run()
    
    tcp     = numpy.load("test_ca_pulse.npz")
    tcp_mod = numpy.load("test_ca_pulse_mod.npz")

    fig, ax = plt.subplots(nrows=4, ncols=1)

    ax[0].plot(tcp_mod['t'], tcp_mod['voltages'])
    ax[0].plot(tcp['t'],     tcp['voltages'],   'r')
    ax[0].set_xlabel("Time [ms]")
    ax[0].set_ylabel("V [mV]")
    ax[0].axis(ymin=-80, ymax=50)

    ax[1].plot(tcp_mod['t'], tcp_mod['ica'])
    ax[1].plot(tcp['t'],     tcp['ica'],   'r')
    ax[1].set_xlabel("Time [ms]")
    ax[1].set_ylabel("ICa [mA/cm2]")
    ax[1].axis(ymin=-1, ymax=0.1)


    ax[2].plot(tcp_mod['t'], tcp_mod['cai'])
    ax[2].plot(tcp['t'], tcp['cai'], 'r')
    ax[2].set_xlabel("Time [ms]")
    ax[2].set_ylabel("Ca [mM]")
    ax[2].axis(ymin=-1E-2, ymax=0.5E-1)

    # ax[1][0].plot(tcp_mod['t'], tcp_mod['cai'] - tcp['cai'])
    # ax[1][0].set_xlabel("Time [ms]")
    # ax[1][0].set_ylabel("Ca [mM]")
    # ax[1][0].axis(ymin=-3E-2, ymax=3E-2)

    ax[3].plot(tcp_mod['t'], tcp_mod['Pi'])
    ax[3].plot(tcp['t'], tcp['Pi'], 'r')
    ax[3].set_xlabel("Time [ms]")
    ax[3].set_ylabel("P [mM]")
    ax[3].axis(ymin=-1E-2, ymax=3E-1)

    # ax[1][1].plot(tcp_mod['t'], tcp_mod['Pi'] - tcp['Pi'])
    # ax[1][1].set_xlabel("Time [ms]")
    # ax[1][1].set_ylabel("P [mM]")
    # ax[1][1].axis(ymin=-3E-2, ymax=3E-2)


    fig.show()
    filename = re.sub('\.', '_', 'compare_ca_pulse-diam%1.1f' % tcp['diam']) + '.pdf'
    fig.savefig('../doc/%s' % filename, format='pdf')

    print('Ca Discrepancy: ' + str(max(abs(tcp_mod['cai'] - tcp['cai']))))
    print('Ca Pc Discrepancy: %2.2f' % (100*max(abs(tcp['cai'] - tcp_mod['cai']))/max(tcp_mod['cai'])))



    print('P Discrepancy: ' + str(max(abs(tcp['Pi'] - tcp_mod['Pi']))))
    print('P Pc Discrepancy: %2.2f' % (100*max(abs(tcp['Pi'] - tcp_mod['Pi']))/max(tcp_mod['Pi'])))

    ## Look just at the times there is an appreciable signal
    #sigmask = (numpy.array(tcp['t']) >= 10) && (numpy.array(tcp['t']) < 22)
    #Pisig = numpy.array(tcp['Pi'][sigmask])
    # Pimodsig = numpy.array(tcp_mod['Pi'][sigmask])
    # print('Mean P disparity in signal: %2.4f' % numpy.mean(Pisig - Pimodsig))
    # print('Mean Ca disparity in signal: %2.4f' % numpy.mean(tcp['Cai'][sigmask] - tcp_mod['Cai'][sigmask]))
