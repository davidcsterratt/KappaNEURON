import matplotlib.pyplot as plt
import matplotlib
import numpy 
from matplotlib.font_manager import FontProperties

dat = numpy.load("simple-psd-pepke-kappa-nmda-comp.npz")
def plot_data(dat, tmax=None, legend=True):
    font = {'family' : 'normal',
            'size'   : 7}
    matplotlib.rc('font', **font)

    times = dat['times'][0]*1E-3
    if (tmax == None): 
        tmax = max(times)
    else:
        tmax = tmax*1E-3

    legendFontP = FontProperties()
    legendFontP.set_size('small')
    fig, (ax1, ax2, ax3, ax4) = plt.subplots(nrows=4, ncols=1, figsize=(2.5, 2.65))
    plt.subplots_adjust(left=0.2, top=0.98, bottom=0.12, right=0.9)

    ax1.plot(times, dat['voltages'][0])
    ax1.set_xlabel("")
    ax1.set_ylabel("V (mV)")
    ax1.axis(ymin=-80, ymax=50)
    ax1.axis(xmin=0, xmax=tmax)
    ax1.xaxis.set_ticklabels([])
    ax1.yaxis.set_ticks([-80, 0])

    ax2.plot(times, dat['ica'][0]*1E3)
    ax2.set_xlabel("")
    ax2.set_ylabel("ICa (uA/cm2)")
    ## ax2.axis(ymin=-1, ymax=0.1)
    ax2.axis(xmin=0, xmax=tmax)
    ax2.xaxis.set_ticklabels([])
    ax2.yaxis.set_ticks([0, -50])

    ax3.plot(times, dat['rec_cai']*1E3)
    ax3.plot(times, dat['rec_CaCBi']*1E3)
    ax3.plot(times, dat['rec_CaCaMNi']*1E3)
    ax3.plot(times, dat['rec_CaCaMCi']*1E3)
    ax3.set_xlabel("")
    ax3.set_ylabel("(uM)")
    plt.axes(ax3)
    if legend:
        plt.legend(('Ca', 'CaCB', 'CaCaMN', 'CaCaMC'), prop=legendFontP)
    ax3.axis(ymin=-1E-3, ymax=0.5E2)
    ## ax3.axis(ymin=-1E-2, ymax=0.5E-1)
    ax3.axis(xmin=0, xmax=tmax)
    ax3.xaxis.set_ticklabels([])
    ax3.yaxis.set_ticks([0, 50])
    
    ax4.plot(times, dat['rec_KCaCaM2Ci']*1E3)
    ax4.plot(times, dat['rec_CaMKIIpi']*1E3)
    # ax4.plot(times, dat['rec_stargazinpi'])
    ax4.set_xlabel("Time (s)")
    ax4.set_ylabel("(uM)")
    plt.axes(ax4)
    if legend:
        plt.legend(('KCaCaM2C','CaMKIIp'), prop=legendFontP)
    ax4.axis(xmin=0, xmax=tmax)
    ax4.axis(ymin=0, ymax=5)
    ax4.yaxis.set_ticks([0, 5])
    #ax4.xaxis.set_ticklabels([0, tmax])

    fig.show() # If the interpreter stops now: close the figure.
    fig.set_size_inches(2.35,2.8)
    # For interactive plotting, see `Part 1` -> `ipython`
    try:
        os.makedirs("figs")
    except Exception as e:
        pass

    fig.savefig("figs/simple-psd-pepke-kappa-nmda-%d.png" % (tmax*1E3), format='png', dpi=600)
    fig.savefig("figs/simple-psd-pepke-kappa-nmda-%d.pdf" % (tmax*1E3))
    return fig

fig1 = plot_data(dat, 600, False)
fig2 = plot_data(dat, 6000)
