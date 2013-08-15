import weakref
import species
import rxdmath
import rxd
import numpy
from generalizedReaction import GeneralizedReaction
from py4j.java_gateway import JavaGateway

gateway = None

class Kappa(GeneralizedReaction):
    def __init__(self, species, kappa_file, regions=None, membrane_flux=False):
        """create a kappa mechanism linked to a species on a given region or set of regions
        if regions is None, then does it on all regions"""
        global gateway
        self._kappa_sims = []
        self._species = weakref.ref(species)
        self._involved_species = [self._species]
        self._kappa_file = kappa_file
        if not hasattr(regions, '__len__'):
            regions = [regions]
        self._regions = regions
        self._trans_membrane = False
        self._update_indices()
        self._membrane_flux = membrane_flux
        if membrane_flux not in (True, False):
            raise Exception('membrane_flux must be either True or False')
        if membrane_flux and regions is None:
            # TODO: rename regions to region?
            raise Exception('if membrane_flux then must specify the (unique) membrane regions')

        if not gateway:
            gateway = JavaGateway()
            print gateway.entry_point

        for i in self._indices[0]:
            print "Creating Kappa Simulation index", i
            kappa_sim = gateway.entry_point.getSpatialKappaSim()
            kappa_sim.loadFile(kappa_file)
            self._kappa_sims.append(kappa_sim)
            ## TODO: Should we check if we are inserting two kappa schemes
            ## in the same place?

        rxd._register_kappa_scheme(self)
    
    def __repr__(self):
        return 'Kappa(%r, kappa_file=%r, regions=%r, membrane_flux=%r)' % (self._species, self._kappa_file, self._regions, self._membrane_flux)
    
    def _update_indices(self):
        # this is called anytime the geometry changes as well as at init
        
        self._indices_dict = {}
        
        # locate the regions containing all species (including the one that changes)
        print self._species()
        if self._species():
            active_regions = [r for r in self._regions if self._species().indices(r)]
        else:
            active_regions = []
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
        
        self._indices = [sum([self._species().indices(r) for r in active_regions], [])]
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

