import os
import weakref
import species
import rxdmath
import rxd
import numpy
import node
from generalizedReaction import GeneralizedReaction
from py4j.java_gateway import JavaGateway

gateway = None

class Kappa(GeneralizedReaction):
    def __init__(self, species, kappa_file, regions=None, membrane_flux=False):
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
        self._update_indices()
        self._membrane_flux = membrane_flux
        if membrane_flux not in (True, False):
            raise Exception('membrane_flux must be either True or False')
        if membrane_flux and regions is None:
            # TODO: rename regions to region?
            raise Exception('if membrane_flux then must specify the (unique) membrane regions')

        rxd._register_kappa_scheme(self)
    
    def __repr__(self):
        return 'Kappa(%r, kappa_file=%r, regions=%r, membrane_flux=%r)' % (self._involved_species, self._kappa_file, self._regions, self._membrane_flux)
    
    def _update_indices(self):
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
            gateway = JavaGateway()
            print gateway.entry_point

        self._kappa_sims = []   # Will this destroy things properly?
        for index in self._indices_dict[self._involved_species[0]()]:
            print "Creating Kappa Simulation in region", r
            kappa_sim = gateway.entry_point.newSpatialKappaSim()
            kappa_sim.loadFile(self._kappa_file)
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
        volumes = node._get_data()[0]
        states = node._get_states()[:]
        for sptr in self._involved_species:
            s = sptr()
            if s:
                for kappa_sim, i in zip(self._kappa_sims, self._indices_dict[s]):
                    nions = round(states[i] \
                                  * rxd._conversion_factor * volumes[i])
                    ## print "Species ", s.name, " conc ", states[i], " nions ", nions
                    kappa_sim.setAgentInitialValue(s.name, nions)
