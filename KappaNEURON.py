import neuron
import neuron.rxd as nr
from neuron.rxd.rxd import *

from scipy.stats import poisson
import math

from py4j.java_gateway import JavaGateway

verbose = False
def report(mess):
    global verbose
    if (verbose):
        print(mess)
_kappa_schemes = []

def _register_kappa_scheme(r):
    # TODO: should we search to make sure that (a weakref to) r hasn't already been added?
    global _kappa_schemes
    _kappa_schemes.append(weakref.ref(r))

def _unregister_kappa_scheme(weakref_r):
    global _kappa_schemes
    _kappa_schemes.remove(weakref_r)

def _kn_setup(): 
    print "Hello World"
    return nr.rxd._setup()
nr.rxd._callbacks[0] = _kn_setup

def _kn_fixed_step_solve(raw_dt):
    import neuron.rxd.rxd as nrr
    print "Step"
    # allow for skipping certain fixed steps
    # warning: this risks numerical errors!
    fixed_step_factor = options.fixed_step_factor
    nrr._fixed_step_count += 1
    if nrr._fixed_step_count % fixed_step_factor: return
    dt = fixed_step_factor * raw_dt
    
    # TODO: this probably shouldn't be here
    if nrr._diffusion_matrix is None and nrr._euler_matrix is None: nrr._setup_matrices()

    states = nrr._node_get_states()[:]

    b = nrr._rxd_reaction(states) - nrr._diffusion_matrix * states
    
    dim = region._sim_dimension
    if dim is None:
        return
    elif dim == 1:
        states[:] += _reaction_matrix_solve(dt, states, _diffusion_matrix_solve(dt, dt * b))

        # clear the zero-volume "nodes"
        states[_zero_volume_indices] = 0

        # TODO: refactor so this isn't in section1d... probably belongs in node
        _section1d_transfer_to_legacy()
    elif dim == 3:
        # the actual advance via implicit euler
        n = len(states)
        m = _scipy_sparse_eye(n, n) - dt * _euler_matrix
        # removed diagonal preconditioner since tests showed no improvement in convergence
        result, info = _scipy_sparse_linalg_bicgstab(m, dt * b)
        assert(info == 0)
        states[:] += result

        for sr in _species_get_all_species().values():
            s = sr()
            if s is not None: s._transfer_to_legacy()

nr.rxd._callbacks[4] = _kn_fixed_step_solve

