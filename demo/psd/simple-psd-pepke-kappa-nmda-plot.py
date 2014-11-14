import matplotlib.pyplot as plt
import numpy 
from matplotlib.font_manager import FontProperties

dat = numpy.load("simple-psd-pepke-kappa-nmda-comp.npz")
def plot_data(dat, tmax=None):
    if (tmax == None): 
        tmax = max(dat['times'][0])
    legendFontP = FontProperties()
    legendFontP.set_size('small')

    fig, (ax1, ax2, ax3, ax4) = plt.subplots(nrows=4, ncols=1, figsize=(2.25*2, 2.65*2))
    plt.subplots_adjust(left=0.2, top=0.98, bottom=0.09)

    ax1.plot(dat['times'][0], dat['voltages'][0])
    ax1.set_xlabel("")
    ax1.set_ylabel("V (mV)")
    ax1.axis(ymin=-80, ymax=50)
    ax1.axis(xmin=0, xmax=tmax)
    ax1.xaxis.set_ticklabels([])
    ax1.yaxis.set_ticks([-80, 0])

    ax2.plot(dat['times'][0], dat['ica'][0])
    ax2.set_xlabel("")
    ax2.set_ylabel("ICa (mA/cm2)")
    ## ax2.axis(ymin=-1, ymax=0.1)
    ax2.axis(xmin=0, xmax=tmax)
    ax2.xaxis.set_ticklabels([])
    ax2.yaxis.set_ticks([0, -0.05])

    ax3.plot(dat['times'][0], dat['rec_cai'])
    ax3.plot(dat['times'][0], dat['rec_CaCBi'])
    ax3.plot(dat['times'][0], dat['rec_CaCaMNi'])
    ax3.plot(dat['times'][0], dat['rec_CaCaMCi'])
    ax3.set_xlabel("")
    ax3.set_ylabel("(mM)")
    plt.axes(ax3)
    plt.legend(('Ca', 'CaCB', 'CaCaMN', 'CaCaMC'), prop=legendFontP)
    ax3.axis(ymin=-1E-3, ymax=0.5E-1)
    ## ax3.axis(ymin=-1E-2, ymax=0.5E-1)
    ax3.axis(xmin=0, xmax=tmax)
    ax3.xaxis.set_ticklabels([])
    ax3.yaxis.set_ticks([0, 0.1])
    
    ax4.plot(dat['times'][0], dat['rec_KCaCaM2Ci'])
    ax4.plot(dat['times'][0], dat['rec_CaMKIIpi'])
    # ax4.plot(dat['times'][0], dat['rec_stargazinpi'])
    ax4.set_xlabel("Time (ms)")
    ax4.set_ylabel("(mM)")
    plt.axes(ax4)
    plt.legend(('KCaCaM2C','CaMKIIp'), prop=legendFontP)
    ax4.axis(xmin=0, xmax=tmax)
    ax4.axis(ymin=0, ymax=0.02)
    ax4.yaxis.set_ticks([0, 0.02])
    ax4.xaxis.set_ticks([0, tmax])

    fig.show() # If the interpreter stops now: close the figure.

    # For interactive plotting, see `Part 1` -> `ipython`
    try:
        os.makedirs("figs")
    except Exception as e:
        pass
        
    fig.savefig("figs/simple-psd-pepke-kappa-nmda-%d.png" % (tmax), format='png', dpi=600)
plot_data(dat, 2000)
plot_data(dat, 20000)
