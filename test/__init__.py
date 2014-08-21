import KappaNEURON
import os
import unittest
import neuron
from neuron import *
from neuron import rxd
import numpy
import matplotlib.pyplot as plt

class TestKappaNEURON(unittest.TestCase):
    ## We can't put this stuff in setUp(), since problems arise if we
    ## try to redefine sections, regions and species. This is because
    ## of class variables in Species.
    ## Set up a kappa simulation of one compartment
    sk = h.Section()
    sk.insert("capulse")
    sk.L=0.2
    sk.diam=0.5
    r = rxd.Region([sk], nrn_region='i')
    ca = rxd.Species(r, name='ca', charge=2, initial=0.0)

    ## Set up a mod simulation of one compartment
    sm = h.Section()
    sm.insert("capulse")
    sm.L=0.2
    sm.diam=0.5
    ## Calcium accumulation via a pump with the pumping turned off
    sm.insert("caPump")
    sm(0.5).k1_caPump = 0 

    def setUp(self):
        self.kappa = KappaNEURON.Kappa([self.ca], self.__module__ + "/caMinimal.ka", self.r, verbose=True)
        KappaNEURON.progress = False
        KappaNEURON.verbose = False

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
        h.dt = h.dt/20
        self.rec_t = h.Vector()
        self.rec_t.record(h._ref_t)
        self.rec_cai = []
        self.rec_v = []
        for sec in h.allsec():
            self.rec_cai.append(h.Vector())
            self.rec_cai[-1].record(sec(0.5)._ref_cai)
            self.rec_v.append(h.Vector())
            self.rec_v[-1].record(sec(0.5)._ref_v)
        init()
        v0 = self.sk(0.5).v
        self.assertEqual(h.t, 0.0)
        self.assertEqual(self.sk(0.5).cai, 0.0)
        run(1.15)
        self.assertAlmostEqual(h.t, 1.15)
        self.assertGreater(self.sk(0.5).cai, 0.0)
        plt.subplots_adjust(left=0.25)
        fig, ax = plt.subplots(nrows=3, ncols=1, figsize=(2.25*2, 2.5*2))
        NA =  6.02214129e23     # Avogadro's number
        caitonum = NA*numpy.pi*(self.sk.diam**2)/4*self.sk.L*1e-18 
        
        i = 0
        for sec in h.allsec():
            ax[0].plot(self.rec_t, self.rec_v[i], color='br'[i])
            ax[1].plot(self.rec_t, self.rec_cai[i], color='br'[i])
            diffv = numpy.diff(self.rec_v[i])
            diffca = caitonum*numpy.diff(self.rec_cai[i])
            ax[2].plot(diffv[1:len(diffv)-1], diffca[0:len(diffv)-2], 'o', color='br'[i])
            fig.show()        
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
            vtocai = cm/(1E-1*2*h.FARADAY*volbyarea)
            print("Eca=%f, t0=%f, t1=%f, gbar=%f, cm=%f, v0=%f, v1=%f" % (eca, t0, t1, gbar, cm, v0, v1))
            print 'Theoretical voltage and Ca difference:'
            Deltav_theo = (eca - v0)*(1 - numpy.exp(-(t1 - t0)*1000*gbar/cm))
            Deltaca_theo = Deltav_theo*vtocai
            print Deltav_theo, Deltaca_theo
            print 'Actual voltage and Ca difference:'
            Deltav = v1 - v0
            Deltaca = sec(0.5).cai
            print Deltav, Deltaca
            self.assertAlmostEqual(Deltav, Deltav_theo, 0)
            self.assertAlmostEqual(Deltaca, Deltaca_theo, 2)
            # self.assertAlmostEqual(v1 - v0, sec(0.5).cai/vtocai, 1)
            i = i + 1

    @unittest.skip("skip for now")
    def test_injectCalcium2(self):
        init()
        self.assertEqual(h.t, 0.0)
        run(3)
        self.assertAlmostEqual(h.t, 3.0)

    def tearDown(self):
        self.kappa = None
        
if __name__ == '__main__':
    unittest.main()
