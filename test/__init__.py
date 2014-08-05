import KappaNEURON
import os
import unittest
from neuron import *
from neuron import rxd

class TestKappaNEURON(unittest.TestCase):
    
    def setUp(self):
        pass

    def test_createSection(self):
        sh = h.Section()
        self.assertIsInstance(sh, nrn.Section)

    def test_injectCalcium(self):
        sh = h.Section()
        r = rxd.Region([sh], nrn_region='i')
        ca = rxd.Species(r, name='ca', charge=2, initial=0.0)
        kappa = KappaNEURON.Kappa([ca], self.__module__ + "/caMinimal.ka", r, verbose=True)
        kappa = []

    def tearDown(self):
        pass
        
if __name__ == '__main__':
    unittest.main()
