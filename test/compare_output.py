import numpy
import matplotlib.pyplot as plt

tcp     = numpy.load("test_ca_pulse.npz")
tcp_mod = numpy.load("test_ca_pulse_mod.npz")

fig, (ax1, ax2, ax3) = plt.subplots(nrows=3, ncols=1)

ax1.plot(tcp['t'], tcp['cai'])
ax1.set_xlabel("Time [ms]")
ax1.set_ylabel("Ca [mM]")
ax1.axis(ymin=-1E-2, ymax=0.5E-1)

ax2.plot(tcp_mod['t'], tcp_mod['cai'])
ax2.set_xlabel("Time [ms]")
ax2.set_ylabel("Ca [mM]")
ax2.axis(ymin=-1E-2, ymax=0.5E-1)

ax3.plot(tcp_mod['t'], tcp_mod['cai'] - tcp['cai'])
ax3.set_xlabel("Time [ms]")
ax3.set_ylabel("Ca [mM]")
ax3.axis(ymin=-3E-2, ymax=3E-2)

fig.show()
