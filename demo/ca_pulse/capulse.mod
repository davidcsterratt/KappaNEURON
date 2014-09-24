TITLE Ca injection pulse
: Sets Ca conductance to gbar between t0 and t1; conductance is 0 otherwise

UNITS {
	  (mA) = (milliamp)
	  (mV) = (millivolt)
    (molar) = (1/liter)
    (mM) = (millimolar)
    FARADAY = (faraday) (coulomb)
    R = (k-mole) (joule/degC)
}

PARAMETER {		:parameters that can be entered when function is called in cell-setup 
	  v             (mV)
	  celsius = 34	(degC)
	  gbar = 1   (mho/cm2)         : initialized conductance
	  cai = 5.e-5   (mM)           : initial internal Ca++ concentration
	  cao = 2       (mM)           : initial external Ca++ concentration
    eca = 140     (mV)           : Ca++ reversal potential
    t0 = 1        (ms)           : Start time of pulse
    t1 = 2        (ms)           : End time of pulse
    fghk = 0                     : Whether GHK equation should be used
}

NEURON {
	  SUFFIX capulse
	  USEION ca READ cai,cao,eca WRITE ica
    RANGE gbar, g, t0, t1, fghk
}

ASSIGNED {                       : parameters needed to solve DE
	  ica   (mA/cm2)
    g  (mho/cm2)
}

INITIAL {                        : initialize the following parameter using rates()
	  g = 0
}

BREAKPOINT {
    g = 0
    if ((t >= t0) && (t <= t1)) {
        g = gbar
    }
    : calcium current induced by this channel
    if (fghk) {
	      ica = g*ghk(v,cai,cao)
    } else {
        ica = g*(v - eca)
    }
}

FUNCTION ghk(v(mV), ci(mM), co(mM)) (mV) {
    LOCAL nu,f
    f = KTF(celsius)/2
    nu = v/f
    ghk=-f*(1. - (ci/co)*exp(nu))*efun(nu)
}

FUNCTION KTF(celsius (degC)) (mV) { : temperature-dependent adjustment factor
    KTF = ((25.(mV)/293.15(degC))*(celsius + 273.15(degC)))
}

FUNCTION efun(z) {
	  if (fabs(z) < 1e-4) {
		    efun = 1 - z/2
	  }else{
		    efun = z/(exp(z) - 1)
	  }
}
