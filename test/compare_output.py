import numpy
import matplotlib.pyplot as plt

tcp     = numpy.load("test_ca_pulse.npz")
tcp_mod = numpy.load("test_ca_pulse_mod.npz")

fig, ax = plt.subplots(nrows=2, ncols=1)

ax[0].plot(tcp_mod['t'], tcp_mod['cai'])
ax[0].set_xlabel("Time [ms]")
ax[0].set_ylabel("Ca [mM]")
ax[0].axis(ymin=-1E-2, ymax=0.5E-1)

ax[0].plot(tcp['t'], tcp['cai'], 'r')
ax[0].set_xlabel("Time [ms]")
ax[0].set_ylabel("Ca [mM]")
ax[0].axis(ymin=-1E-2, ymax=0.5E-1)

# ax[1][0].plot(tcp_mod['t'], tcp_mod['cai'] - tcp['cai'])
# ax[1][0].set_xlabel("Time [ms]")
# ax[1][0].set_ylabel("Ca [mM]")
# ax[1][0].axis(ymin=-3E-2, ymax=3E-2)

ax[1].plot(tcp_mod['t'], tcp_mod['Pi'])
ax[1].set_xlabel("Time [ms]")
ax[1].set_ylabel("P [mM]")
ax[1].axis(ymin=-1E-2, ymax=3E-1)

ax[1].plot(tcp['t'], tcp['Pi'], 'r')
ax[1].set_xlabel("Time [ms]")
ax[1].set_ylabel("P [mM]")
ax[1].axis(ymin=-1E-2, ymax=3E-1)

# ax[1][1].plot(tcp_mod['t'], tcp_mod['Pi'] - tcp['Pi'])
# ax[1][1].set_xlabel("Time [ms]")
# ax[1][1].set_ylabel("P [mM]")
# ax[1][1].axis(ymin=-3E-2, ymax=3E-2)


fig.show()

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

