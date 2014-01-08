import neuron
import neuron.rxd as nr

def _kn_setup(): 
    print "Hello World"
    return nr.rxd._setup()
nr.rxd._callbacks[0] = _kn_setup


