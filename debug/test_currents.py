from neuron import *
from neuron import rxd
import neuron.rxd.rxd as nrr
import random
import matplotlib.pyplot as plt

def _kn_currents(rhs):
    nrr._currents(rhs)
    cur = -random.random()/10
    ## This line alters ica, and cai but on its own does not seem to
    ## affect the voltage
    nrr._curr_ptrs[0][0] += cur
    ## This line is necessary to change the voltage
    rhs[1] += -cur

nrr._callbacks[2] = _kn_currents

sh = h.Section()
sh.insert("pas")
r = rxd.Region([sh], nrn_region='i')
ca = rxd.Species(r, name='ca', charge=2, initial=0.0)

rec_t = h.Vector()
rec_t.record(h._ref_t)
rec_v = h.Vector()
rec_v.record(sh(0.5)._ref_v)
rec_ica = h.Vector()
rec_ica.record(sh(0.5)._ref_ica)
rec_cai = h.Vector()
rec_cai.record(sh(0.5)._ref_cai)

init()
run(15)

fig, ax = plt.subplots(nrows=3, ncols=1, figsize=(2.25*2, 2.5*2))
plt.subplots_adjust(left=0.2, top=0.95, bottom=0.1)

ax[0].plot(rec_t, rec_v)
ax[0].set_xlabel("")
ax[0].set_ylabel("V (mV)")
ax[0].axis(ymin=-80, ymax=50)

ax[1].plot(rec_t, rec_ica)
ax[1].set_xlabel("")
ax[1].set_ylabel("ICa (mA/cm2)")
ax[1].axis(ymin=-0.1, ymax=0.1)

ax[2].plot(rec_t, rec_cai)
ax[2].set_xlabel("")
ax[2].set_ylabel("[Ca] (mM)")
ax[2].axis(ymin=0, ymax=0.001)

fig.show()

