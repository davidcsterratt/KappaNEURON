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
    print "Setup"
    return nr.rxd._setup()
    global _kappa_schemes
    # update Kappa schemes
    for kptr in _kappa_schemes:
        k = kptr()
        if k is not None: k.re_init()

nr.rxd._callbacks[0] = _kn_setup

def _kn_fixed_step_solve(raw_dt):
    import neuron.rxd.rxd as nrr
    global _kappa_schemes
    
    report("---------------------------------------------------------------------------")
    report("FIXED STEP SOLVE. NEURON time %f" % h.t)
    report("states")

    # allow for skipping certain fixed steps
    # warning: this risks numerical errors!
    fixed_step_factor = options.fixed_step_factor
    nrr._fixed_step_count += 1
    if nrr._fixed_step_count % fixed_step_factor: return
    dt = fixed_step_factor * raw_dt
    
    # TODO: this probably shouldn't be here
    if nrr._diffusion_matrix is None and nrr._euler_matrix is None: nrr._setup_matrices()

    states = nrr._node_get_states()[:]
    report(states)

    report("flux b")
    ## DCS: This gets fluxes (from ica, ik etc) and computes changes
    ## due to reactions
    b = nrr._rxd_reaction(states) - nrr._diffusion_matrix * states
    report(b)
    
    dim = region._sim_dimension
    if dim is None:
        return
    elif dim == 1:
        states[:] += _reaction_matrix_solve(dt, states, _diffusion_matrix_solve(dt, dt * b))

        ## Go through each kappa scheme. The region belonging to each
        ## kappa scheme should not overlap with any other kappa scheme's
        ## region.
        volumes = node._get_data()[0]
        for kptr in _kappa_schemes:
            k = kptr()

            ## Now we want add any fluxes to the kappa sims and update the
            ## quantities seen in NEURON.

            ## There is one kappa_sim for each active region in the kappa
            ## scheme.

            report("\nRUN 0.5 KAPPA STEP")
            for kappa_sim in k._kappa_sims:
                ## kappa_sim.runByTime2(h.t - dt/2)      # Second argument is "time per
                kappa_sim.runByTime3(dt/2)      # Second argument is "time per

            ## This should work for multiple species working, but has only
            ## been tested for ca
            report("\nADDING FLUXES TO KAPPA")
            for  sptr in k._involved_species:
                s = sptr()
                name = s.name
                report("ION: %s" % (name))
                for kappa_sim, i in zip(k._kappa_sims, k._indices_dict[s]):
                    ## Number of ions
                    ## Flux b has units of mM/ms
                    ## Volumes has units of um3
                    ## _conversion factor has units of molecules mM^-1 um^-3
                    mu = dt * b[i] * _conversion_factor * volumes[i]
                    nions = 0.0
                    if mu!=0:
                        nions = math.copysign(1, mu)*poisson.rvs(abs(mu))
                    report("index %d; volume: %f ; flux %f ; # of ions: %s" % (i, volumes[i], b[i], nions))
                    kappa_sim.addAgent(name, nions)
                    t_kappa = kappa_sim.getTime()
                    discrepancy = h.t - t_kappa
                    report('Kappa Time %f; NEURON time %f; Discrepancy %f' % (t_kappa, h.t, discrepancy))


            report("\nRUN 0.5 KAPPA STEP")  
            for kappa_sim in k._kappa_sims:
                ## kappa_sim.runByTime2(h.t)      # Second argument is "time per
                kappa_sim.runByTime3(dt/2)      # Second argument is "time per
                t_kappa = kappa_sim.getTime()
                discrepancy = h.t - t_kappa
                report('Kappa Time %f; NEURON time %f; Discrepancy %f' % (t_kappa, h.t, discrepancy))
                ## if (abs(discrepancy) > 1e-3):
                ##    raise NameError('NEURON time (%f) does not match Kappa time (%f). Discrepancy = %f ' % (h.t, t_kappa, h.t - t_kappa))

            ## Update states
            for  sptr in k._involved_species:
                s = sptr()
                name = s.name
                for kappa_sim, i in zip(k._kappa_sims, k._indices_dict[s]):
                    states[i] = kappa_sim.getObservation(name) \
                        /(_conversion_factor * volumes[i])

        report("Updated states")
        report(states)

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

