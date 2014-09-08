import KappaNEURON
import os
import unittest
import neuron
from neuron import *
from neuron import rxd
import numpy as np
import matplotlib.pyplot as plt
import re

class TestCaAccumulation(unittest.TestCase):
    ## Whether to plot
    plot = True

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
    
    ## Time of pulse
    t0 = 1.0
    t1 = 1.1
    tstop = t1

    ## Conversions
    NA =  6.02214129e23     # Avogadro's number

    gbar = 0.001
    cm = sk(0.5).cm

    k1 = 0                      # Pump parameter
    k2 = False                  # Optional parameter
    P = False              # Pump species
    P0 = 0 

    tol = 0.01
    def assertEqualWithinTol(self, a, b, tol=None):
        if tol == None:
            tol = self.tol
        self.assertAlmostEqual(a, b, delta=tol*a)
        

    def setUp(self):
        self.caitonum = self.NA*np.pi*(self.sk.diam**2)/4*self.sk.L*1e-18 
        h.dt = h.dt/4

    def get_mode(self, sec):
        ## Determine if section contains mod pump or kappa pump
        mode = 'kappa'
        for mech in sec(0.5):
            if re.search('caPump', mech.name()):
                mode = 'mod'
        return(mode)

    def injectCalcium(self, ghk=0, mechanism='caPump1'):
        
        ## Insert calcium pump into mod section
        self.sm.insert(mechanism)

        ## Insert calcium pump into kappa section
        if mechanism == 'caPump1':
            self.kappa = KappaNEURON.Kappa([self.ca], self.__module__ + "/" + mechanism + ".ka", self.r, verbose=True)
        if mechanism == 'caPump2':
            self.P  = rxd.Species(self.r, name='P', charge=0, initial=self.P0)
            self.kappa = KappaNEURON.Kappa([self.ca, self.P], self.__module__ + "/" + mechanism + ".ka", self.r, verbose=True)
            self.kappa.setVariable('vol', self.sk.L*(self.sk.diam**2)/4*np.pi)
            self.kappa.setVariable('k2', self.k2)
            setattr(self.sm(0.5), 'k2_' + mechanism, self.k2)

        ## Set variables
        self.kappa.setVariable('k1', self.k1)
        setattr(self.sm(0.5), 'k1_' + mechanism, self.k1)

        ## Set up recordings
        self.rec_t = h.Vector()
        self.rec_t.record(h._ref_t)
        self.rec_cai = []
        self.rec_v = []
        for sec in h.allsec():
            self.rec_cai.append(h.Vector())
            self.rec_cai[-1].record(sec(0.5)._ref_cai)
            self.rec_v.append(h.Vector())
            self.rec_v[-1].record(sec(0.5)._ref_v)

        # self.sm(0.5).P0_caPump2 = self.P0
        if mechanism == 'caPump2':
            self.rec_Pi = []
            for sec in h.allsec():
                self.rec_Pi.append(h.Vector())
                print self.get_mode(sec)
                if self.get_mode(sec) == 'mod':
                    self.rec_Pi[-1].record(sec(0.5)._ref_P_caPump2)
                else:
                    self.rec_Pi[-1].record(sec(0.5)._ref_Pi)
        
        KappaNEURON.progress = False
        KappaNEURON.verbose = False

        self.assertIsInstance(self.sk, nrn.Section)
        self.assertEqual(self.ca.initial, 0.0)
        for sec in h.allsec():
            ## This forces eca to be a constant, rather than being
            ## computed from Nernst equation at every time step
            ## See http://www.neuron.yale.edu/neuron/static/new_doc/modelspec/programmatic/ions.html?highlight=ion_style#ion_style
            if (ghk == 1):
                sec.fghk_capulse = 1
            else:
                sec.push(); h('ion_style("ca_ion", 3, 1, 0, 0, 1)') ; h.pop_section()
                
            for seg in sec:
                seg.t0_capulse = self.t0
                seg.t1_capulse = self.t1

                seg.gbar_capulse = self.gbar

        ## Initialise simulation
        init()
        eca0 = self.sk(0.5).eca
        self.v0 = self.sk(0.5).v
        self.assertEqual(h.t, 0.0)
        self.assertEqual(self.sk(0.5).cai, 0.0)
        ## Run
        run(self.tstop)
        self.assertAlmostEqual(h.t, self.tstop)
        self.assertGreater(self.sk(0.5).cai, 0.0)
        if (ghk == 0):
            self.assertEqual(eca0, self.sk(0.5).eca)


    def get_diffv_diffca(self, i):
        times = np.array(self.rec_t)
        stim_inds = np.where((times > self.t0) & (times < self.t1))
        ## Check that during the stimulus, every voltage increment
        ## is proportional to the calcium increment in the
        ## *preceeding* timestep.  At the end of the pulse and the
        ## beginning of the pulse this is not true, because the voltage increment 
        diffv = np.diff(np.array(self.rec_v[i])[stim_inds])
        diffca = np.diff(np.array(self.rec_cai[i])[stim_inds])
        return(diffv, diffca)

    def do_plot(self):
        plt.subplots_adjust(left=0.25)
        nrow = 3
        if (not self.P):
            nrow = 4
        fig, ax = plt.subplots(nrows=4, ncols=1, figsize=(2.25*2, 2.5*2))
        i = 0
        for sec in h.allsec():
            (diffv, diffca) = self.get_diffv_diffca(i)
            ax[0].plot(self.rec_t, self.rec_v[i], color='br'[i])
            (Deltav_theo, Deltaca_theo) = self.get_Deltav_Deltaca_theo(sec, np.array(self.rec_t))
            ax[0].plot(self.rec_t, Deltav_theo + self.v0, color='g')
            ax[1].plot(self.rec_t, self.rec_cai[i], color='br'[i])
            ax[2].plot(diffv[1:len(diffv)-1], self.caitonum*diffca[0:len(diffv)-2], 'o', color='br'[i])
            if (self.P0 > 0):
                ax[3].plot(self.rec_t, self.rec_Pi[i], color='br'[i])
            fig.show()        
            i = i + 1

    def get_Deltav_Deltaca_theo(self, sec, t, verbose=False):
        eca = sec(0.5).eca
        volbyarea = sec.diam/4
        vtocai = self.cm/(1E-1*2*h.FARADAY*volbyarea)
        ## Print some variables
        v1 = sec(0.5).v
        if verbose:
            print("cca=%f, t0=%f, t1=%f, gbar=%f, cm=%f, v0=%f, v1=%f" % (sec(0.5).eca, self.t0, self.t1, self.gbar, self.cm, self.v0, v1))
        if verbose:
            print 'Theoretical voltage and Ca difference:'
        tau = self.cm/(1000.0*self.gbar + self.k1*self.cm)
        vinf = (1000.0*self.gbar*eca + self.cm*self.k1*self.v0)/(1000.0*self.gbar + self.cm*self.k1)
        Deltav_theo = (vinf - self.v0)*(1-np.exp(-(t - self.t0)/tau))
        if isinstance(Deltav_theo, numpy.ndarray):
            Deltav_theo[np.where(t < self.t0)] = 0.0
        Deltaca_theo = Deltav_theo*vtocai
        if verbose:
            print "eca= %f; tau=%f; vinf=%f; Deltav_theo=%f; Deltaca_theo=%f" % (eca, tau, vinf, Deltav_theo, Deltaca_theo)
        return(Deltav_theo, Deltaca_theo)

    def get_Deltav_Deltaca(self, sec, i):
        v1 = sec(0.5).v
        print 'Actual voltage and Ca difference:'
        Deltav = v1 - self.v0
        Deltaca = sec(0.5).cai - self.rec_cai[i][0]
        print Deltav, Deltaca
        return(Deltav, Deltaca)

    def get_stats(self):
        Deltav       = {}
        Deltaca      = {}
        Deltav_theo  = {}
        Deltaca_theo = {}
        volbyarea    = {}
        vtocai       = {}
        diffv        = {}
        diffca       = {}
        
        i = 0
        for sec in h.allsec():
            ## Print some variables
            v1 = sec(0.5).v
            print("Eca=%f, t0=%f, t1=%f, gbar=%f, cm=%f, v0=%f, v1=%f" % (sec(0.5).eca, self.t0, self.t1, self.gbar, self.cm, self.v0, v1))

            ## Determine if section contains mod pump or kappa pump
            mode = self.get_mode(sec)
            print mode
            (Deltav_theo[mode],  Deltaca_theo[mode])  = self.get_Deltav_Deltaca_theo(sec, self.tstop)
            (Deltav[mode],       Deltaca[mode])       = self.get_Deltav_Deltaca(sec, i)

            volbyarea[mode] = sec.diam/4
            vtocai[mode] = self.cm/(1E-1*2*h.FARADAY*volbyarea[mode])

            (diffv[mode], diffca[mode]) = self.get_diffv_diffca(i)
            i = i+1

        return(Deltav, Deltaca, Deltav_theo, Deltaca_theo, volbyarea, vtocai, diffv, diffca)


    def test_injectCalcium(self):
        self.tstop = self.t1 + h.dt
        self.injectCalcium(ghk=0)
        self.do_plot()
        ## Needed for tests involving theoretical voltage and calcium to work
        self.tstop = self.t1
        ## Get data from both sections
        Deltav, Deltaca, Deltav_theo, Deltaca_theo, volbyarea, vtocai, diffv, diffca = self.get_stats()

        ## Check theory and simulation match, if using
        ## deterministic ('mod') simulation
        self.assertAlmostEqual(Deltav['mod'], Deltav_theo['mod'], 0)
        self.assertAlmostEqual(Deltaca['mod'], Deltaca_theo['mod'], 2)

        ## Calcium and voltage should be in sync
        for mode in ['mod', 'kappa']:
            self.assertAlmostEqual(Deltav[mode], Deltaca[mode]/vtocai[mode], 0)

        ## All calcium ion increments should be integers
        self.assertAlmostEqual(max(self.caitonum*diffca['kappa'] - np.round(self.caitonum*diffca['kappa'])), 0, 2)
        ## Calcium ion increments should be equal to voltage increments
        self.assertAlmostEqual(max(abs(vtocai['kappa']*diffv['kappa'][1:len(diffv['kappa'])-1] - diffca['kappa'][0:len(diffv['kappa'])-2])), 0, 2)


    def test_injectCalciumGHK(self):
        self.injectCalcium(ghk=1)
        self.do_plot()
        ## Get data from both sections
        Deltav, Deltaca, Deltav_theo, Deltaca_theo, volbyarea, vtocai, diffv, diffca = self.get_stats()

        ## Calcium and voltage should be in sync
        for mode in ['mod', 'kappa']:
            self.assertAlmostEqual(Deltav[mode], Deltaca[mode]/vtocai[mode], 0)

        ## All calcium ion increments should be integers
        self.assertAlmostEqual(max(self.caitonum*diffca['kappa'] - np.round(self.caitonum*diffca['kappa'])), 0, 2)
        ## Calcium ion increments should be equal to voltage increments
        self.assertAlmostEqual(max(abs(vtocai['kappa']*diffv['kappa'][1:len(diffv['kappa'])-1] - diffca['kappa'][0:len(diffv['kappa'])-2])), 0, 2)

    def test_injectCalciumPump(self):
        self.t1 = 2
        self.tstop = 2
        self.k1 = 1
        self.injectCalcium(ghk=0)
        self.do_plot()
        ## Get data from both sections
        Deltav, Deltaca, Deltav_theo, Deltaca_theo, volbyarea, vtocai, diffv, diffca = self.get_stats()

        ## Check theory and simulation match, if using
        ## deterministic ('mod') simulation
        self.assertAlmostEqual(Deltav['mod'],  Deltav_theo['mod'],  0)
        self.assertAlmostEqual(Deltaca['mod'], Deltaca_theo['mod'], 2)

        ## Calcium and voltage should be in sync
        for mode in ['mod', 'kappa']:
            self.assertEqualWithinTol(Deltav[mode], Deltaca[mode]/vtocai[mode])

        ## Check that kappa and deterministic simulations agree to
        ## within 10%
        tol = 0.1
        self.assertLess(abs((Deltav['kappa'] - Deltav['mod'])/(Deltav['mod'] - self.v0)), 0.1)
        self.assertLess(abs((Deltaca['kappa'] - Deltaca['mod'])/Deltaca['mod']), 0.1)
        

    def test_injectCalciumPumpGHK(self):
        self.t1 = 2
        self.tstop = 2
        self.k1 = 1
        self.injectCalcium(ghk=1)
        self.do_plot()
        ## Get data from both sections
        Deltav, Deltaca, Deltav_theo, Deltaca_theo, volbyarea, vtocai, diffv, diffca = self.get_stats()

        ## Calcium and voltage should be in sync
        for mode in ['mod', 'kappa']:
            self.assertAlmostEqual(Deltav[mode], Deltaca[mode]/vtocai[mode], 0)

        ## Check that kappa and deterministic simulations agree to
        ## within 10%
        tol = 0.1
        self.assertLess(abs((Deltav['kappa'] - Deltav['mod'])/(Deltav['mod'] - self.v0)), 0.1)
        self.assertLess(abs((Deltaca['kappa'] - Deltaca['mod'])/Deltaca['mod']), 0.1)

        ## All calcium ion increments should be integers
        self.assertAlmostEqual(max(self.caitonum*diffca['kappa'] - np.round(self.caitonum*diffca['kappa'])), 0, 2)
        ## Calcium ion increments should be equal to voltage increments
        self.assertAlmostEqual(max(abs(vtocai['kappa']*diffv['kappa'][1:len(diffv['kappa'])-1] - diffca['kappa'][0:len(diffv['kappa'])-2])), 0, 2)


    def test_injectCalciumPump2(self):
        self.t1 = 2.0
        self.tstop = 3.0
        self.k1 = 1
        self.k2 = 0
        self.P0 = 0.20 
        self.injectCalcium(ghk=0, mechanism='caPump2')
        self.do_plot()

        ## Run through both sections
        times = np.array(self.rec_t)
        i = 0
        for sec in h.allsec():
            mode = self.get_mode(sec)
            print mode
            v = np.array(self.rec_v[i])
            self.assertAlmostEqual(v[np.where(np.isclose(times, self.t1 + 0.1))],
                                   v[np.where(np.isclose(times, self.tstop))])
            cai = np.array(self.rec_cai[i])
            self.assertGreater(cai[np.where(np.isclose(times, self.t1 + 0.1))],
                               cai[np.where(np.isclose(times, self.tstop))])


    def tearDown(self):
        self.kappa = None

testSuite = unittest.TestSuite()
#testSuite.addTest(TestCaAccumulation('test_injectCalcium'))
testSuite.addTest(TestCaAccumulation('test_injectCalciumGHK'))
        
if __name__ == '__main__':
    unittest.main()
