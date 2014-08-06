import KappaNEURON
import os
import unittest
import neuron
from neuron import *
from neuron import rxd

class TestKappaNEURON(unittest.TestCase):
    def setUp(self):
        pass

    def test_createSection(self):
        print 'test_createSection'
        sh = h.Section()
        self.assertIsInstance(sh, nrn.Section)
        sh.insert("capulse")

    def test_injectCalcium(self):
        print 'test_injectCalcium'
        sh = h.Section()
        sh.insert("capulse")
        r = rxd.Region([sh], nrn_region='i')
        ca = rxd.Species(r, name='ca', charge=2, initial=0.0)
        kappa = KappaNEURON.Kappa([ca], self.__module__ + "/caMinimal.ka", r, verbose=True)
        init()
        run(30)
        # print h.t, h.v
        kappa = []

    def tearDown(self):
        pass
        
if __name__ == '__main__':
    unittest.main()
