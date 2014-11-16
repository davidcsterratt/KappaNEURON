## Spine head into which Ca flows via a very simple channel and is
## pumped out.  Voltage clamp ensures almost constant Ca flow.
from neuron import *
from subprocess import call
import re
import matplotlib.pyplot as plt
import matplotlib
import numpy

vinit = -70

def make_spine_head(L=0.1, diam=0.2, gbar=0.05):
    # Spine Head
    sh = h.Section()
    sh.insert("pas")                # Passive channel
    sh.insert("capulse")            # Code to give Ca pulse
    sh.L = L
    sh.diam = diam

    ## This setting of parameters gives a calcium influx and pump
    ## activation that is more-or-less scale-independent
    sh.gbar_capulse = gbar # *sh.diam
    sh.fghk_capulse = 1         # Use GHK
    sh.t0_capulse = 5
    sh.t1_capulse = 10
    
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
    global vinit
    init()
    h.finitialize(vinit)
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

def plot_records(tcp_mod, tcp):
    font = {'family' : 'normal',
            'size'   : 7}
    matplotlib.rc('font', **font)
    # figsize=(2.25, 3)
    fig, ax = plt.subplots(nrows=4, ncols=1, figsize=(2.35,2.8))
    plt.subplots_adjust(left=0.25, top=0.95, bottom=0.15,right=0.95)

    ax[0].plot(tcp['t'],     tcp['voltages'],     'r')
    ax[0].plot(tcp_mod['t'], tcp_mod['voltages'], 'b')
    ax[0].set_xlabel("")
    ax[0].set_ylabel("V (mV)")
    ax[0].axis(ymin=-100, ymax=0)
    ax[0].xaxis.set_ticklabels([])
    ax[0].yaxis.set_ticks([-100, 0])

    icamin = min(numpy.concatenate([tcp_mod['ica'], tcp['ica']]))
    icamax = max(numpy.concatenate([tcp_mod['ica'], tcp['ica']]))
    ax[1].plot(tcp['t'],     tcp['ica'],     'r')
    ax[1].plot(tcp_mod['t'], tcp_mod['ica'], 'b')
    ax[1].set_xlabel("")
    ax[1].set_ylabel("ICa (mA/cm2)")
    ax[1].axis(ymin=icamin*1.1,
               ymax=icamax*1.1)
    ax[1].xaxis.set_ticklabels([])
    ax[1].yaxis.set_ticks([0.1, 0, -0.2])

    ## Convert [Ca] in mM to uM
    caimax = max(numpy.concatenate([tcp_mod['cai'], tcp['cai']])*1000.0)
    ax[2].plot(tcp['t'],     tcp['cai']*1000.0,     'r')
    ax[2].plot(tcp_mod['t'], tcp_mod['cai']*1000.0, 'b')
    ax[2].set_xlabel("")
    ax[2].set_ylabel("[Ca] (uM)")
    ax[2].axis(ymin=-0.1*caimax, ymax=caimax*1.1)
    ax[2].xaxis.set_ticklabels([])
    ##ax[2].yaxis.set_ticks([0, 0.020])

    Pimax = max(numpy.concatenate([tcp_mod['Pi'], tcp['Pi']]))
    ax[3].plot(tcp['t'],     tcp['Pi'],     'r')
    ax[3].plot(tcp_mod['t'], tcp_mod['Pi'], 'b')
    ax[3].set_xlabel("Time (ms)")
    ax[3].set_ylabel("[P] (mM)")
    ax[3].axis(ymin=0, ymax=Pimax*1.1)
    ax[3].yaxis.set_ticks([0, 0.20])

    fig.show()
    fig.set_size_inches(2.35,2.6)
    return fig, ax

def compare_traces(diam=0.2, gbar=0.001, 
                   gamma2=0.1, P0=0.2, vclamp=False):
    import test_ca_pulse_mod
    test_ca_pulse_mod.run(diam=diam, 
                          gbar=gbar,
                          gamma2=gamma2,
                          P0=P0,
                          vclamp=vclamp)
    import test_ca_pulse
    test_ca_pulse.run(diam=diam, 
                      gbar=gbar,
                      gamma2=gamma2,
                      P0=P0,
                      vclamp=vclamp)
    
    tcp     = numpy.load("test_ca_pulse.npz")
    tcp_mod = numpy.load("test_ca_pulse_mod.npz")

    fig, ax = plot_records(tcp_mod, tcp)

    filename = re.sub('\.', '_', 'compare_ca_pulse-diam:%1.1f-gbar:%1.3f-gamma2:%1.3f-P0:%1.3f' % (diam, gbar, gamma2, P0)) + '.pdf'
    fig.savefig('%s' % filename, format='pdf')

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

def animate_traces():
    tcp     = numpy.load("test_ca_pulse.npz")
    tcp_mod = numpy.load("test_ca_pulse_mod.npz")

    fig, ax = plot_records(tcp_mod, tcp)
    Tmax = int(numpy.floor(max(tcp['t'])));
    Vlim = [-80, 50]
    Ilim = [-0.02, 0.001]
    cailim = [-1E-2, 0.2E-1]
    Pilim = [-1E-2, 1E-1]
    for T in numpy.nditer(numpy.arange(0, Tmax, 0.5)):
        inds = numpy.array(tcp_mod['t']) <= T
        ax[0].cla()
        ax[0].set_xlim([0, Tmax])
        ax[0].set_ylim(Vlim)
        ax[0].plot(tcp_mod['t'][inds], tcp_mod['voltages'][inds])
        ax[0].set_ylabel("V [mV]")

        ax[1].cla()
        ax[1].set_xlim([0, Tmax])
        ax[1].set_ylim(Ilim)
        ax[1].plot(tcp_mod['t'][inds], tcp_mod['ica'][inds])
        ax[1].set_ylabel("ICa [mA/cm2]")

        ax[2].cla()
        ax[2].set_xlim([0, Tmax])
        ax[2].set_ylim(cailim)
        ax[2].plot(tcp_mod['t'][inds], tcp_mod['cai'][inds])
        ax[2].set_ylabel("Ca [mM]")

        ax[3].cla()
        ax[3].set_xlim([0, Tmax])
        ax[3].set_ylim(Pilim)
        ax[3].plot(tcp_mod['t'][inds], tcp_mod['Pi'][inds])
        ax[3].set_ylabel("P [mM]")

        print(T)
        filename = 'animation/test_ca_pulse%04d.png' % (T*10)
        print(filename)
        fig.savefig(filename)

    for T in numpy.nditer(numpy.arange(0, Tmax, 0.5)):
        inds = numpy.array(tcp_mod['t']) <= T
        ax[0].cla()
        ax[0].set_xlim([0, Tmax])
        ax[0].set_ylim(Vlim)
        ax[0].plot(tcp_mod['t'], tcp_mod['voltages'])
        ax[0].plot(tcp['t'][inds], tcp['voltages'][inds], 'r')
        ax[0].set_ylabel("V [mV]")

        ax[1].cla()
        ax[1].set_xlim([0, Tmax])
        ax[1].set_ylim(Ilim)
        ax[1].plot(tcp_mod['t'], tcp_mod['ica'])
        ax[1].plot(tcp['t'][inds], tcp['ica'][inds], 'r')
        ax[1].set_ylabel("ICa [mA/cm2]")

        ax[2].cla()
        ax[2].set_xlim([0, Tmax])
        ax[2].set_ylim(cailim)
        ax[2].plot(tcp_mod['t'], tcp_mod['cai'])
        ax[2].plot(tcp['t'][inds], tcp['cai'][inds], 'r')
        ax[2].set_ylabel("Ca [mM]")

        ax[3].cla()
        ax[3].set_xlim([0, Tmax])
        ax[3].set_ylim(Pilim)
        ax[3].plot(tcp_mod['t'], tcp_mod['Pi'])
        ax[3].plot(tcp['t'][inds], tcp['Pi'][inds], 'r')
        ax[3].set_ylabel("P [mM]")

        fig.savefig('animation/test_ca_pulse%04d.png' % ((T + Tmax)*10))

    os.system("mencoder 'mf://animation/*.png' -mf type=png:fps=10 -ovc lavc -lavcopts vcodec=wmv2 -oac copy -o animation.mpg")

        


