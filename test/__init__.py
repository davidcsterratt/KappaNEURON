import KappaNEURON
import os
import unittest
import neuron
from neuron import *
from neuron import rxd

class TestKappaNEURON(unittest.TestCase):
    def setUp(self):
        self.sh = h.Section()
        self.sh.insert("capulse")
        self.sh.L=2
        self.sh.diam=1
        self.r = rxd.Region([self.sh], nrn_region='i')
        if 'ca' not in rxd.species._get_all_species():
            self.ca = rxd.Species(self.r, name='ca', charge=2, initial=0.0)
        else:
            self.ca = rxd.species._get_all_species()['ca']()
        self.kappa = KappaNEURON.Kappa([self.ca], self.__module__ + "/caMinimal.ka", self.r, verbose=True)
        KappaNEURON.progress = False
        self.rec_cai = h.Vector()
        self.rec_cai.record(self.sh(0.5)._ref_cai)

    def test_createSection(self):
        self.assertIsInstance(self.sh, nrn.Section)
        init()
        self.assertEqual(h.t, 0.0)
        self.assertEqual(self.sh(0.5).cai, 0.0)

    def test_injectCalcium(self):
        init()
        self.assertEqual(h.t, 0.0)
        # self.assertEqual(self.sh(0.5).cai, 0.0)
        run(30)

    def test_injectCalcium2(self):
        init()
        run(30)

    def tearDown(self):
        self.kappa = None
        
if __name__ == '__main__':
    unittest.main()
