import KappaNEURON
import os
import unittest
import neuron
from neuron import *
from neuron import rxd
import numpy

class TestKappaNEURON(unittest.TestCase):
    ## We can't put this stuff in setUp(), since problems arise if we
    ## try to redefine sections, regions and species. This is because
    ## of class variables in Species.
    ## Set up a kappa simulation of one compartment
    sk = h.Section()
    sk.insert("capulse")
    sk.L=2
    sk.diam=1
    r = rxd.Region([sk], nrn_region='i')
    ca = rxd.Species(r, name='ca', charge=2, initial=0.0)

    ## Set up a mod simulation of one compartment
    sm = h.Section()
    sm.insert("capulse")
    sm.L=2
    sm.diam=1
    ## Calcium accumulation via a pump with the pumping turned off
    sm.insert("caPump")
    sm(0.5).k1_caPump = 0 

    def setUp(self):
        self.kappa = KappaNEURON.Kappa([self.ca], self.__module__ + "/caMinimal.ka", self.r, verbose=True)
        KappaNEURON.progress = False
        KappaNEURON.verbose = False
        self.rec_cai = h.Vector()
        self.rec_cai.record(self.sk(0.5)._ref_cai)

    def test_injectCalcium(self):
        self.assertIsInstance(self.sk, nrn.Section)
        self.assertEqual(self.ca.initial, 0.0)
        for sec in h.allsec():
            ## This forces eca to be a constant, rather than being
            ## computed from Nernst equation at every time step
            ## See http://www.neuron.yale.edu/neuron/static/new_doc/modelspec/programmatic/ions.html?highlight=ion_style#ion_style
            sec.push(); h('ion_style("ca_ion", 3, 1, 0, 0, 1)') ; h.pop_section()
            for seg in sec:
                seg.t1_capulse = 1.1
                seg.gbar_capulse = 0.001
        # h.dt = h.dt/20
        init()
        v0 = self.sk(0.5).v
        self.assertEqual(h.t, 0.0)
        self.assertEqual(self.sk(0.5).cai, 0.0)
        run(3)
        self.assertAlmostEqual(h.t, 3.0)
        self.assertGreater(self.sk(0.5).cai, 0.0)
        for sec in h.allsec():
            print sec.name()
            for mech in sec(0.5):
                print mech.name()
            v1 = sec(0.5).v
            volbyarea = self.sk.diam/4
            eca = sec(0.5).eca
            t0 = sec(0.5).t0_capulse
            t1 = sec(0.5).t1_capulse
            gbar = sec(0.5).gbar_capulse
            cm = sec(0.5).cm
            print("Eca=%f, t0=%f, t1=%f, gbar=%f, cm=%f, v0=%f" % (eca, t0, t1, gbar, cm, v0))
            print 'Theoretical voltage and Ca difference:'
            Deltav_theo = (eca - v0)*(1 - numpy.exp(-(t1 - t0)*1000*gbar/cm))
            Deltaca_theo = Deltav_theo/(1E-1*2*h.FARADAY*volbyarea)*cm
            print Deltav_theo, Deltaca_theo
            print 'Actual voltage and Ca difference:'
            Deltav = v1 - v0
            Deltaca = sec(0.5).cai
            print Deltav, Deltaca
            self.assertAlmostEqual(Deltav, Deltav_theo, 0)
            self.assertAlmostEqual(Deltaca, Deltaca_theo, 2)
            # self.assertAlmostEqual(v1 - v0, 1E-1*2*h.FARADAY*volbyarea*sec(0.5).cai/sec(0.5).cm)

    @unittest.skip("skip for now")
    def test_injectCalcium2(self):
        init()
        self.assertEqual(h.t, 0.0)
        run(30)
        self.assertAlmostEqual(h.t, 30.0)

    def tearDown(self):
        self.kappa = None
        
if __name__ == '__main__':
    unittest.main()
