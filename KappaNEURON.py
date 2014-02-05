import neuron
import neuron.rxd as nr
import neuron.rxd.rxd as nrr

from scipy.stats import poisson
import math

import SpatialKappa
from py4j.protocol import * 

import re

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

def _kn_init(): 
    print "_kn_init()"
    nrr._init()
    global _kappa_schemes
    # update Kappa schemes
    for kptr in _kappa_schemes:
        k = kptr()
        if k is not None: k.re_init()

#
# register the initialization handler and the advance handler
#
nrr._fih = neuron.h.FInitializeHandler(_kn_init)

## Override the NEURON nonvint _fixed_step_solve callback   
def _kn_fixed_step_solve(raw_dt):
    global _kappa_schemes
    
    report("---------------------------------------------------------------------------")
    report("FIXED STEP SOLVE. NEURON time %f" % nrr.h.t)
    report("states")

    # allow for skipping certain fixed steps
    # warning: this risks numerical errors!
    fixed_step_factor = nrr.options.fixed_step_factor
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

    ## DCS FIXME: This is different from the old rxd.py file - need check what
    ## the difference is
    b = nrr._rxd_reaction(states) - nrr._diffusion_matrix * states
    report(b)
    
    dim = nrr.region._sim_dimension
    if dim is None:
        return
    elif dim == 1:
        states[:] += nrr._reaction_matrix_solve(dt, states, nrr._diffusion_matrix_solve(dt, dt * b))

        ## Go through each kappa scheme. The region belonging to each
        ## kappa scheme should not overlap with any other kappa scheme's
        ## region.
        volumes = nrr.node._get_data()[0]
        for kptr in _kappa_schemes:
            k = kptr()

            ## Now we want add any fluxes to the kappa sims and update the
            ## quantities seen in NEURON.

            ## There is one kappa_sim for each active region in the kappa
            ## scheme.

            report("\nRUN 0.5 KAPPA STEP")
            for kappa_sim in k._kappa_sims:
                kappa_sim.runForTime(dt/2)      # Second argument is "time per

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
                    mu = dt * b[i] * nrr._conversion_factor * volumes[i]
                    nions = 0.0
                    if mu!=0:
                        nions = numpy.sign(mu)*poisson.rvs(abs(mu))
                    report("index %d; volume: %f ; flux %f ; # of ions: %s" % (i, volumes[i], b[i], nions))
                    kappa_sim.addAgent(name, nions)
                    t_kappa = kappa_sim.getTime()
                    discrepancy = nrr.h.t - t_kappa
                    report('Kappa Time %f; NEURON time %f; Discrepancy %f' % (t_kappa, nrr.h.t, discrepancy))


            report("\nRUN 0.5 KAPPA STEP")  
            for kappa_sim in k._kappa_sims:
                kappa_sim.runForTime(dt/2)      # Second argument is "time per
                t_kappa = kappa_sim.getTime()
                discrepancy = nrr.h.t - t_kappa
                report('Kappa Time %f; NEURON time %f; Discrepancy %f' % (t_kappa, nrr.h.t, discrepancy))
                ## if (abs(discrepancy) > 1e-3):
                ##    raise NameError('NEURON time (%f) does not match Kappa time (%f). Discrepancy = %f ' % (h.t, t_kappa, h.t - t_kappa))

            ## Update states
            for  sptr in k._involved_species:
                s = sptr()
                name = s.name
                for kappa_sim, i in zip(k._kappa_sims, k._indices_dict[s]):
                    states[i] = kappa_sim.getObservation(name) \
                        /(nrr._conversion_factor * volumes[i])

        report("Updated states")
        report(states)

        # clear the zero-volume "nodes"
        states[nrr._zero_volume_indices] = 0

        # TODO: refactor so this isn't in section1d... probably belongs in node
        nrr._section1d_transfer_to_legacy()
    elif dim == 3:
        # the actual advance via implicit euler
        n = len(states)
        m = _scipy_sparse_eye(n, n) - dt * _euler_matrix
        # removed diagonal preconditioner since tests showed no improvement in convergence
        result, info = _scipy_sparse_linalg_bicgstab(m, dt * b)
        assert(info == 0)
        states[:] += result

        for sr in nrr._species_get_all_species().values():
            s = sr()
            if s is not None: s._transfer_to_legacy()

nrr._callbacks[4] = _kn_fixed_step_solve

import os
import weakref
import neuron.rxd.species
import neuron.rxd.rxdmath
import neuron.rxd
import numpy
import neuron.rxd.node
from neuron.rxd.generalizedReaction import GeneralizedReaction

gateway = None

class Kappa(GeneralizedReaction):
    def __init__(self, species, kappa_file, regions=None, membrane_flux=False, time_units='ms', verbose=False):
        """create a kappa mechanism linked to a species on a given region or set of regions
        if regions is None, then does it on all regions"""
        global gateway
        self._kappa_sims = []
        self._species = []
        for s in species:
            self._species.append(weakref.ref(s))
        ## self._species = weakref.ref(species)
        self._involved_species = self._species
        self._kappa_file = os.path.join(os.getcwd(), kappa_file)
        print(os.getcwd())
        if not hasattr(regions, '__len__'):
            regions = [regions]
        self._regions = regions
        self._active_regions = []
        self._trans_membrane = False
        self._membrane_flux = membrane_flux
        self._time_units = 'ms'
        self._time_units = time_units
        self._verbose = verbose
        if membrane_flux not in (True, False):
            raise Exception('membrane_flux must be either True or False')
        if membrane_flux and regions is None:
            # TODO: rename regions to region?
            raise Exception('if membrane_flux then must specify the (unique) membrane regions')
        self._update_indices()
        print('Registering kappa scheme')
        _register_kappa_scheme(self)
        print _kappa_schemes
        self._weakref = weakref.ref(self) # Seems to be needed for the destructor
    
    def __repr__(self):
        return 'Kappa(%r, kappa_file=%r, regions=%r, membrane_flux=%r)' % (self._involved_species, self._kappa_file, self._regions, self._membrane_flux)
    
    def __del__(self):
        ## A similar idiom to rxd._register_kappa_scheme() doesn't seem to work
        _unregister_kappa_scheme(self._weakref)
        for kappa_sim in self._kappa_sims:
            del(kappa_sim)

    def _update_indices(self):
        print '_update_indices'
        global gateway

        # this is called anytime the geometry changes as well as at init
        
        self._indices_dict = {}
        
        # locate the regions containing all species (including the one
        # that channges)

        active_regions = self._regions
        for sptr in self._involved_species:
            s = sptr()
            if s:
                for r in self._regions:
                    if r in active_regions and not s.indices(r):
                        del active_regions[active_regions.index(r)]
            else:
                active_regions = []
        
        # store the indices
        for sptr in self._involved_species:
            s = sptr()
            self._indices_dict[s] = sum([s.indices(r) for r in active_regions], [])
        ## Check that each species has the same number of elements
        if (len(set([len(self._indices_dict[s()]) for s in self._involved_species])) != 1):
            raise Exception('Different numbers of indices for various species') 
        self._active_regions = active_regions

        ## Create the kappa simulations
        if not gateway:
            gateway = SpatialKappa.SpatialKappa()
            print gateway

        self._kappa_sims = []   # Will this destroy things properly?
        for index in self._indices_dict[self._involved_species[0]()]:
            print "Creating Kappa Simulation in region", r
            kappa_sim = gateway.kappa_sim(self._time_units, verbose)
            try:
                kappa_sim.loadFile(self._kappa_file)
            except Py4JJavaError as e:
                java_err = re.sub(r'java.lang.IllegalStateException: ', r'', str(e.java_exception))
                errstr = 'Error in kappa file %s: %s' % (self._kappa_file, java_err)
                raise RuntimeError(errstr)
                
            self._kappa_sims.append(kappa_sim)
            ## TODO: Should we check if we are inserting two kappa schemes
            ## in the same place?

        self._mult = [1]

    def _do_memb_scales(self):
        # TODO: does anyone still call this?
        # TODO: update self._memb_scales (this is just a dummy value to make things run)
        self._memb_scales = 1


    
    def _get_memb_flux(self, states):
        if self._membrane_flux:
            raise Exception('membrane flux due to rxd.Rate objects not yet supported')
            # TODO: refactor the inside of _evaluate so can construct args in a separate function and just get self._rate() result
            rates = self._evaluate(states)[2]
            return self._memb_scales * rates
        else:
            return []

    def setVariable(self, variable, value):
        for kappa_sim in self._kappa_sims:
            kappa_sim.setVariable(float(value), variable)

    ## This is perhaps an abuse of this function, but it is called at
    ## init() time
    def re_init(self):
        print "Kappa.re_init()"
        volumes = nrr.node._get_data()[0]
        states = nrr.node._get_states()[:]
        print(states)
        for sptr in self._involved_species:
            s = sptr()
            if s:
                for kappa_sim, i in zip(self._kappa_sims, self._indices_dict[s]):
                    nions = round(states[i] \
                                  * nrr._conversion_factor * volumes[i])
                    ## print "Species ", s.name, " conc ", states[i], " nions ", nions
                    try:
                        kappa_sim.setAgentInitialValue(s.name, nions)
                    except:
                        print('Error setting initial value of agent %s to %d' % (s.name, nions))
                        raise


    def run_free(self, t_run):
        # Run free of neuron
        for kptr in _kappa_schemes:
            k = kptr()
            for kappa_sim in k._kappa_sims:
                t_kappa = kappa_sim.getTime()
                kappa_sim.runForTime(t_kappa + t_run)
