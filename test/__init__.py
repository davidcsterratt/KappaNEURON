import KappaNEURON
import os
import unittest
import neuron
from neuron import *
from neuron import rxd

class TestKappaNEURON(unittest.TestCase):
    ## We can't put this stuff in setUp(), since problems arise if we
    ## try to redefine sections, regions and species. This is because
    ## of class variables in Species.
    sh = h.Section()
    sh.insert("capulse")
    sh.L=2
    sh.diam=1
    r = rxd.Region([sh], nrn_region='i')
    ca = rxd.Species(r, name='ca', charge=2, initial=0.0)

    def setUp(self):
        self.kappa = KappaNEURON.Kappa([self.ca], self.__module__ + "/caMinimal.ka", self.r, verbose=True)
        KappaNEURON.progress = False
        KappaNEURON.verbose = False
        self.rec_cai = h.Vector()
        self.rec_cai.record(self.sh(0.5)._ref_cai)

    def test_injectCalcium(self):
        self.assertIsInstance(self.sh, nrn.Section)
        self.assertEqual(self.ca.initial, 0.0)
        init()
        self.assertEqual(h.t, 0.0)
        self.assertEqual(self.sh(0.5).cai, 0.0)
        run(15)
        self.assertAlmostEqual(h.t, 15.0)
        self.assertGreater(self.sh(0.5).cai, 0.0)

    def test_injectCalcium2(self):
        init()
        self.assertEqual(h.t, 0.0)
        run(30)
        self.assertAlmostEqual(h.t, 30.0)

    def tearDown(self):
        self.kappa = None
        
if __name__ == '__main__':
    unittest.main()
