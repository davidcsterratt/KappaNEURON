import neuron
from neuron import h
import neuron.rxd as nr
import neuron.rxd.rxd as nrr
import neuron.rxd.species
import neuron.rxd.rxdmath
import neuron.rxd.node
from neuron.rxd.generalizedReaction import GeneralizedReaction, molecules_per_mM_um3
from neuron.rxd.multiCompartmentReaction import MultiCompartmentReaction
import weakref
import random
import itertools

import SpatialKappa
from py4j.protocol import * 

from scipy.stats import poisson
import numpy
import re
import os, sys
import warnings

FARADAY = h.FARADAY

verbose = False
def report(mess):
    global verbose
    if (verbose):
        print(mess)
_kappa_schemes = []

progress = True

def _register_kappa_scheme(r):
    # TODO: should we search to make sure that (a weakref to) r hasn't already been added?
    global _kappa_schemes
    _kappa_schemes.append(weakref.ref(r))

def _unregister_kappa_scheme(weakref_r):
    global _kappa_schemes
    _kappa_schemes.remove(weakref_r)

def _kn_init(): 
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

mode = 'lumped_influx' # 'continuous_influx'
mode = 'continuous_influx'

def _kn_fixed_step_solve(raw_dt):
    if (mode == 'lumped_influx'):
        _kn_fixed_step_solve_lumped_influx(raw_dt)
    else:
        _kn_fixed_step_solve_continuous_influx(raw_dt)

## Override the NEURON nonvint _fixed_step_solve callback   
def _kn_fixed_step_solve_lumped_influx(raw_dt):
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
                kappa_sim.runForTime(dt/2, False)      # Second argument is "time per

            ## This should work for multiple species working, but has only
            ## been tested for ca
            report("\nADDING FLUXES TO KAPPA")
            for  s in k._membrane_species:
                name = s.name
                report("ION: %s" % (name))
                for kappa_sim, i in zip(k._kappa_sims, k._indices_dict[s]):
                    ## Number of ions
                    ## Flux b has units of mM/ms
                    ## Volumes has units of um3
                    ## _conversion factor has units of molecules mM^-1 um^-3
                    mu = dt * b[i] * molecules_per_mM_um3 * volumes[i]
                    nions = 0
                    if mu!=0:
                        nions = numpy.sign(mu)*poisson.rvs(abs(mu))
                    report("index %d; volume: %f ; flux %f ; # ions: %s; mu: %f\n" % (i, volumes[i], b[i], nions, mu))
                    try:
                        kappa_sim.addAgent(name, nions)
                    except Py4JJavaError as e:
                        java_err = re.sub(r'java.lang.IllegalStateException: ', r'', str(e.java_exception))
                        errstr = 'Problem adding agents, probably trying to add too many in one step; try reducing number of agents by reducing surface area or current density:\n%s' % (java_err)
                        raise RuntimeError(errstr)

                    t_kappa = kappa_sim.getTime()
                    discrepancy = nrr.h.t - t_kappa
                    report('Kappa Time %f; NEURON time %f; Discrepancy %f' % (t_kappa, nrr.h.t, discrepancy))


            report("\nRUN 0.5 KAPPA STEP")  
            for kappa_sim in k._kappa_sims:
                kappa_sim.runForTime(dt/2, False)      # Second argument is "time per
                t_kappa = kappa_sim.getTime()
                discrepancy = nrr.h.t - t_kappa + dt/2
                report('Kappa Time %f; NEURON time %f; Discrepancy %f' % (t_kappa, nrr.h.t, discrepancy))
                ## This code is commented out because it doesn't work
                ## if run_free has been used; this makes NEURON and
                ## SpatialKappa time go out of sync
                ## 
                ## if (abs(discrepancy) > 1e-3):
                ##     raise NameError('NEURON time (%f) does not match Kappa time (%f). Discrepancy = %f ' % (nrr.h.t + dt/2, t_kappa, discrepancy))

            ## Update states
            for  sptr in k._involved_species:
                s = sptr()
                name = s.name
                for kappa_sim, i in zip(k._kappa_sims, k._indices_dict[s]):
                    states[i] = kappa_sim.getVariable(name) \
                        /(molecules_per_mM_um3 * volumes[i])

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
    
    t = nrr.h.t + dt
    if progress:
        sys.stdout.write("\rTime = %12.5f/%5.5f [%3.3f%%]" % (t, neuron.h.tstop, t/neuron.h.tstop*100))
        if (abs(t - neuron.h.tstop) < 1E-6):
            sys.stdout.write("\n")
        sys.stdout.flush()

def _run_kappa_continuous(states, b, dt):
    #############################################################################
    ## 1. Pass all relevant continous variables to the rule-based simulator
    ##
    ## Relevant variables might be
    ## * Calcium current (for deterministic channels)
    ## * Membrane potential (for stochastic channels controlled by Kappa model)
    #############################################################################

    ## Go through each kappa scheme. The region belonging to each
    ## kappa scheme should not overlap with any other kappa scheme's
    ## region.
    volumes = nrr.node._get_data()[0]
    for kptr in _kappa_schemes:
        k = kptr()
        report("\nPASSING FLUXES TO KAPPA")
        for s in k._membrane_species:
            report("ION: %s" % (s.name))
            for kappa_sim, i in zip(k._kappa_sims, k._indices_dict[s]):
                ## Number of ions
                ## Flux b has units of mM/ms
                ## Volumes has units of um3
                ## _conversion factor has units of molecules mM^-1 um^-3
                flux = b[i] * molecules_per_mM_um3 * volumes[i]
                kappa_sim.setTransitionRate('Create %s' % (s.name), flux)
                report("Setting %s flux[%d] to b[%d]*NA*vol[%d] = %f*%f*%f = %f" % (s.name, i, i, i,  b[i], molecules_per_mM_um3, volumes[i], flux))

            report("PASSING MEMBRANE POTENTIAL TO KAPPA")
            ## TODO: pass membrane potential to kappa

        #############################################################################
        ## 2. Run the rule-based simulator from t to t + dt
        #############################################################################
        #############################################################################
        ## 3. Compute the net change Delta Stot in the number of each
        ## bridging species S and convert back into a current.
        #############################################################################
        #############################################################################
        ## 4. Set the corresponding elements of the flux to the
        ## currents computed in step 3
        #############################################################################

        for kptr in _kappa_schemes:
            k = kptr()

            ## Recording total starting value of each species
            Stot0 = {}
            for f in k._kappa_fluxes:
                for sptr in f._sources:
                    s = sptr()._species()
                    Stot0[s.name] = {}
                    for kappa_sim, i in zip(k._kappa_sims, k._indices_dict[s]):
                        Stot0[s.name][i] = kappa_sim.getVariable('Total %s' % (s.name))
                        report("Stot0[%s][%d] = %f" % (s.name, i, Stot0[s.name][i]))

            report("RUN 1 KAPPA STEP")  
            for kappa_sim in k._kappa_sims:
                kappa_sim.runForTime(dt, False)     
                t_kappa = kappa_sim.getTime()
                report("kappa time now %f" % (t_kappa))

            ## Recording total ending value of each membrane species
            for f in k._kappa_fluxes:
                f._memb_flux = []
                for sptr in f._sources:
                    s = sptr()._species()
                    for kappa_sim, i in zip(k._kappa_sims, k._indices_dict[s]):
                        ## For ions, compute the current
                        Stot1 = kappa_sim.getVariable('Total %s' % (s.name))
                        report("Stot1[%s][%d] = %f" % (s.name, i, Stot1))
                        DeltaStot = Stot1 - Stot0[s.name][i]
                        bnew = DeltaStot/(dt*molecules_per_mM_um3*volumes[i])
                        f._memb_flux.append(-(bnew - b[i]))
                        report("Species %s: DeltaStot=%d, bnew=%f, b=%f, _memb_flux=%f" % (s.name, DeltaStot, bnew, b[i], f._memb_flux[-1]))
                        b[i] = bnew
        
        report("States before update")
        report(states)
        
        #############################################################################
        ## 5. Update the continous variables according to the update step
        #############################################################################
        states[:] += nrr._reaction_matrix_solve(dt, states, nrr._diffusion_matrix_solve(dt, dt * b))
        report("States after continuous update")
        report(states)

        #############################################################################
        ## 6. Voltage step overrides states, possibly making them negative so put back actual states
        #############################################################################
        for kptr in _kappa_schemes:
            k = kptr()
            ## Recording total ending value of each species
            for sptr in k._involved_species:
                s = sptr()
                for kappa_sim, i in zip(k._kappa_sims, k._indices_dict[s]):
                    ## Update concentration
                    states[i] = kappa_sim.getVariable(s.name) \
                                /(molecules_per_mM_um3 * volumes[i])
        report("States after kappa update")
        report(states)

    return states


## Override the NEURON nonvint _fixed_step_solve callback   
def _kn_fixed_step_solve_continuous_influx(raw_dt):
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
        states = _run_kappa_continuous(states, b, dt)

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
    
    t = nrr.h.t + dt
    sys.stdout.write("\rTime = %12.5f/%5.5f [%3.3f%%]" % (t, neuron.h.tstop, t/neuron.h.tstop*100))
    if (abs(t - neuron.h.tstop) < 1E-6):
        sys.stdout.write("\n")
    sys.stdout.flush()


nrr._callbacks[4] = _kn_fixed_step_solve

gateway = None

def setSeed(seed):
    global _kappa_schemes
    print _kappa_schemes
    k = _kappa_schemes[0]
    k()._kappa_sims[0].setSeed(seed)

    # _kappa_sims[0].setSeed(seed)

class Kappa(GeneralizedReaction):
    def __init__(self, *args, **kwargs):
        """Specify a Kappa model spanning the membrane to be added to the system.
        
        Use this for, for example, pumps and channels, or interactions between
        species living in a volume (e.g. the cytosol) and species on a
        membrane (e.g. the plasma membrane).

        Keyword arguments:

        membrane_species -- List of rxd.Species defining which species
        in the Kappa model cross the cell membrane

        species -- List of rxd.Species defining species to observe
        inside the Kappa model.

        kappa_file -- Name of a Kappa file defining rules. The file
        should contain agents with the same names as all the
        membrane_species and observables with names corresponding to
        names of the rxd.Species in the species argument.

        regions -- List of rxd.Regions in which the Kappa mechanism
        should be inserted.

        membrane_flux -- Boolean indicating if the reaction should
        produce a current across the plasma membrane that should
        affect the membrane potential.
        
        .. seealso::
        
            :class:`neuron.rxd.multiCompartmentReaction`

        """
        
        # additional keyword arguments
        membrane_species = kwargs.get('membrane_species', [])
        species = kwargs.get('species', [])
        kappa_file = kwargs.get('kappa_file')
        regions = kwargs.get('regions', None)
        membrane_flux = kwargs.get('membrane_flux', True)
        time_units = kwargs.get('time_units', 'ms')
        scale_by_area = kwargs.get('scale_by_area', True)

        global gateway
        self._kappa_sims = []
        self._species = []
        for s in membrane_species + species:
            self._species.append(weakref.ref(s))
            if s.initial is None:
                s.initial = 0
                warnings.warn('Initial concentration of %s not specified; setting to zero' % (s.name), UserWarning)
        ## This is the species that crosses the membrane
        self._membrane_species = membrane_species
        ## self._species = weakref.ref(species)
        self._involved_species = self._species
        self._kappa_file = os.path.join(os.getcwd(), kappa_file)
        if not hasattr(regions, '__len__'):
            regions = [regions]
        self._regions = regions
        self._active_regions = []
        self._scale_by_area = scale_by_area
        self._trans_membrane = False
        self._membrane_flux = False
        self._time_units = time_units
        if membrane_flux not in (True, False):
            raise Exception('membrane_flux must be either True or False')
        if membrane_flux and regions is None:
            raise Exception('if membrane_flux then must specify the (unique) membrane regions')
        ## Create KappaFlux objects
        self._kappa_fluxes = []
        for s in self._membrane_species:
            self._kappa_fluxes.append(KappaFlux(membrane_species=[s], 
                                                regions=self._regions,
                                                membrane_flux=membrane_flux,
                                                kappa_parent=self))
        self._sources = []
        self._dests = []
        self._update_indices()
        self._create_kappa_sims()
        report('Registering kappa scheme')
        _register_kappa_scheme(self)
        nrr._register_reaction(self)
        report(_kappa_schemes)
        self._weakref = weakref.ref(self) # Seems to be needed for the destructor
    
    def __repr__(self):
        return 'Kappa(%r, kappa_file=%r, regions=%r, membrane_flux=%r)' % (self._involved_species, self._kappa_file, self._regions, self._membrane_flux)
    
    def __del__(self):
        global gateway, _kappa_schemes
        ## A similar idiom to rxd._register_kappa_scheme() doesn't seem to work
        _unregister_kappa_scheme(self._weakref)
        for kappa_sim in self._kappa_sims:
            del(kappa_sim)
        ## Needed to ensure cleanup and no exit errors in python2.7
        if (len(_kappa_schemes) == 0):
            gateway = None

    def _create_kappa_sims(self):
        """Create the kappa simulations."""
        global gateway

        ## Create the kappa simulations
        if not gateway:
            gateway = SpatialKappa.SpatialKappa()

        self._kappa_sims = []   # Will this destroy things properly?
        for index in self._indices_dict[self._involved_species[0]()]:
            report("Creating Kappa Simulation in index %d" % (index))
            kappa_sim = gateway.kappa_sim(self._time_units, verbose)
            try:
                kappa_sim.loadFile(self._kappa_file)
            except Py4JJavaError as e:
                java_err = re.sub(r'java.lang.IllegalStateException: ', r'', str(e.java_exception))
                errstr = 'Error in kappa file %s: %s' % (self._kappa_file, java_err)
                raise RuntimeError(errstr)
                
            ## Set up transitions to create membrane species in Kappa
            ## simulation and measure the total amount of the agent
            ## corresponding to the membrane species
            if (mode == 'continuous_influx'):
                for s in self._membrane_species:
                    ## Get description of agent
                    agent = kappa_sim.getAgentMap(s.name)
                    link_names = agent[s.name].keys()
                    if (len(link_names) > 1):
                        errstr = 'Error in kappa file %s: Agent %s has more than one site' % (self._kappa_file, s.name)
                        raise RuntimeError()
                    
                    link_name = link_names[0]
                    
                    ## Add transition to create 
                    kappa_sim.addTransition('Create %s' % (s.name), {}, agent, 0.0)
                    
                    ## Add variable to measure total species
                    kappa_sim.addVariableMap('Total %s' % (s.name), {s.name: {link_name: {'l': '?'}}})
                    ## Add observation variable
                    kappa_sim.addVariableMap('%s' % (s.name), {s.name: {}})

            self._kappa_sims.append(kappa_sim)
            ## TODO: Should we check if we are inserting two kappa schemes
            ## in the same place?
        self._mult = [1]

    def setVariable(self, variable, value):
        """Sets a variable in the Kappa simuluations.

        Keyword arguments:

        variable -- String corresponding to %var described in kappa_file.

        value -- Float to set the variable to.

        """
        for kappa_sim in self._kappa_sims:
            kappa_sim.setVariable(float(value), variable)

    def run_free(self, t_run):
        """Run Kappa simulations free of NEURON

        During run_free() invocations, there is no passing or
        receiving of fluxes between NEURON and Kappa.

        Keyword arguments:

        t_run -- Time in millseconds for which to run.

        """
        for kptr in _kappa_schemes:
            k = kptr()
            for kappa_sim in k._kappa_sims:
                kappa_sim.runForTime(float(t_run), True)

    ## Overridden functions
    def re_init(self):
        """Sets the initial concentration/number of kappa variables.  

        This is perhaps an abuse of this function, but it is called at
        init() time.

        """
        volumes = nrr.node._get_data()[0]
        states = nrr.node._get_states()[:]
        for sptr in self._involved_species:
            s = sptr()
            if s:
                for kappa_sim, i in zip(self._kappa_sims, self._indices_dict[s]):
                    nions = round(states[i] \
                                  * molecules_per_mM_um3 * volumes[i])
                    ## print "Species ", s.name, " conc ", states[i], " nions ", nions
                    try:
                        kappa_sim.getVariable(s.name)
                    except:
                        raise NameError('There is no observable or variable in %s called %s; add a line like this:\n%%obs: \'%s\' <complex definition> ' % (self._kappa_file, s.name, s.name))
                    if kappa_sim.isAgent(s.name):
                        try:
                            kappa_sim.setAgentInitialValue(s.name, nions)
                        except Py4JJavaError as e:
                            raise NameError('Error setting initial value of agent %s to %d\n%s' % (s.name, nions,  str(e.java_exception)))

    def _evaluate(self, states):
        """This does nothing in the KappaNEURON class"""
        return ([], [], [])

    def _jacobian_entries(self, states, multiply=1, dx=1.e-10):
        """This does nothing in the KappaNEURON class"""
        return ([], [], [])        


class KappaFlux(MultiCompartmentReaction):
    def __init__(self, *args, **kwargs):
        """Specify a KappaFlux spanning the membrane to be added to the system.
        
        Use this for, for example, pumps and channels, or interactions between
        species living in a volume (e.g. the cytosol) and species on a
        membrane (e.g. the plasma membrane).

        Keyword arguments:

        membrane_species -- List of rxd.Species defining which species
        in the Kappa model cross the cell membrane

        species -- List of rxd.Species defining species to observe
        inside the Kappa model.

        regions -- List of rxd.Regions in which the Kappa mechanism
        should be inserted.

        membrane_flux -- Boolean indicating if the reaction should
        produce a current across the plasma membrane that should
        affect the membrane potential.
        
        .. seealso::
        
            :class:`neuron.rxd.multiCompartmentReaction`

        """
        
        # additional keyword arguments
        membrane_species = kwargs.get('membrane_species', [])
        regions = kwargs.get('regions', None)
        kappa_parent = kwargs.get('kappa_parent', None)
        membrane_flux = kwargs.get('membrane_flux', True)
        scale_by_area = kwargs.get('scale_by_area', True)

        self._species = []
        for s in membrane_species:
            self._species.append(weakref.ref(s))
            if s.initial is None:
                s.initial = 0
                warnings.warn('Initial concentration of %s not specified; setting to zero' % (s.name), UserWarning)
        ## This is the species that crosses the membrane
        self._membrane_species = membrane_species
        ## self._species = weakref.ref(species)
        self._involved_species = self._species
        self._kappa_parent = kappa_parent
        if not hasattr(regions, '__len__'):
            regions = [regions]
        self._regions = regions
        self._active_regions = []
        self._scale_by_area = scale_by_area
        self._trans_membrane = False
        self._membrane_flux = membrane_flux
        if membrane_flux not in (True, False):
            raise Exception('membrane_flux must be either True or False')
        if membrane_flux and regions is None:
            raise Exception('if membrane_flux then must specify the (unique) membrane regions')
        self._memb_flux = None
        ## Set up the sources for _get_memb_flux(). In
        ## multicompartmentReaction.py some of this is done in
        ## _update_rates()
        self._lhs_items = []
        self._sources = []
        for s in self._membrane_species:
            i = s[self._regions[0]]
            self._lhs_items.append(i)
            w = weakref.ref(i)
            self._sources += [w]
        self._dests = []
        self._update_indices()
        nrr._register_reaction(self)
        self._weakref = weakref.ref(self) # Seems to be needed for the destructor

    def __repr__(self):
        return 'KappaFlux(%r, kappa_parent=%r, regions=%r, membrane_flux=%r)' % (self._involved_species, self._kappa_parent, self._regions, self._membrane_flux)

    def _evaluate(self, states):
        """This does nothing in the KappaNEURON class"""
        return ([], [], [])

    def _jacobian_entries(self, states, multiply=1, dx=1.e-10):
        """This does nothing in the KappaNEURON class"""
        return ([], [], [])        
    
    def _get_memb_flux(self, states):
        """Returns the flux across the membrane due to univalent ion in mA/cm^2

        In KappaNEURON, this flux is determined by the net change in
        number of ions during the preceding time step, which is
        calculated in _kn_fixed_step_solve().
        """
        if self._membrane_flux:
            ## _memb_flux has been set in _kn_fixed_step_solve(), unless it's
            ## the first time step, in which case we need to create it.
            if self._memb_flux is None:
                len_memb_flux = 0
                for sptr in self._sources:
                    s = sptr()._species()
                    len_memb_flux += len(self._indices_dict[s])
                    self._memb_flux = nrr._numpy_zeros(len_memb_flux)

            ## TODO: Use the full volumes and surface_area vectors
            volumes, surface_area, diffs = nrr.node._get_data()
            ## This has units of mA/cm2
            return -self._memb_scales*self._memb_flux*molecules_per_mM_um3
        else:
            return []

    def _do_memb_scales(self, cur_map): 
        """Set up self._memb_scales and cur_map."""
        if not self._scale_by_area:
            areas = numpy.ones(len(areas))
        else:
            volumes = numpy.concatenate([list(self._regions[0]._geometry.volumes1d(sec) for sec in self._regions[0].secs)])
        neuron_areas = []
        for sec in self._regions[0].secs:
            neuron_areas += [h.area((i + 0.5) / sec.nseg, sec=sec) for i in xrange(sec.nseg)]
        neuron_areas = numpy.array(neuron_areas)
        # area_ratios is usually a vector of 1s
        area_ratios = volumes / neuron_areas
        # still needs to be multiplied by the valence of each molecule
        self._memb_scales = -area_ratios * FARADAY / (10000 * molecules_per_mM_um3)
        #print area_ratios
        #print self._memb_scales
        #import sys
        #sys.exit()
        
        # since self._memb_scales is only used to compute currents as seen by the rest of NEURON,
        # we only use NEURON's areas 
        #self._memb_scales = volume * molecules_per_mM_um3 / areas
        
        
        if self._membrane_flux:
            # TODO: don't assume/require always inside/outside on one side...
            #       if no nrn_region specified, then (make so that) no contribution
            #       to membrane flux
            source_regions = [s()._region()._nrn_region for s in self._sources]
            dest_regions = [d()._region()._nrn_region for d in self._dests]
            
            if 'i' in source_regions and 'o' not in source_regions and 'i' not in dest_regions:
                inside = -1 #'source'
            elif 'o' in source_regions and 'i' not in source_regions and 'o' not in dest_regions:
                inside = 1 # 'dest'
            elif 'i' in dest_regions and 'o' not in dest_regions and 'i' not in source_regions:
                inside = 1 # 'dest'
            elif 'o' in dest_regions and 'i' not in dest_regions and 'o' not in source_regions:
                inside = -1 # 'source'
            else:
                raise RxDException('unable to identify which side of reaction is inside (hope to remove the need for this')
        
        # dereference the species to get the true species if it's actually a SpeciesOnRegion
        sources = [s()._species() for s in self._sources]
        dests = [d()._species() for d in self._dests]
        if self._membrane_flux:
            if any(s in dests for s in sources) or any(d in sources for d in dests):
                # TODO: remove this limitation
                raise RxDException('current fluxes do not yet support same species on both sides of reaction')
        
        # TODO: make so don't need multiplicity (just do in one pass)
        # TODO: this needs changed when I switch to allowing multiple sides on the left/right (e.g. simplified Na/K exchanger)
        self._cur_charges = tuple([-inside * s.charge for s in sources if s.name is not None] + [inside * s.charge for s in dests if s.name is not None])
        self._net_charges = sum(self._cur_charges)
        
        self._cur_ptrs = []
        self._cur_mapped = []
        
        for sec in self._regions[0].secs:
            for i in xrange(sec.nseg):
                local_ptrs = []
                local_mapped = []
                for sp in itertools.chain(self._sources, self._dests):
                    spname = sp()._species().name
                    if spname is not None:
                        name = '_ref_i%s' % (spname)
                        seg = sec((i + 0.5) / sec.nseg)
                        local_ptrs.append(seg.__getattribute__(name))
                        uberlocal_map = [None, None]
                        if spname + 'i' in cur_map:
                            uberlocal_map[0] = cur_map[spname + 'i'][seg]
#                        if spname + 'o' in cur_map:
#                            uberlocal_map[1] = cur_map[spname + 'o'][seg]
                        local_mapped.append(uberlocal_map)
                self._cur_ptrs.append(tuple(local_ptrs))
                self._cur_mapped.append(tuple(local_mapped))

