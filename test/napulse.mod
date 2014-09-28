TITLE Sodium injection pulse
: Sets Glu current to inabar between t0 and t1; conductance is 0 otherwise

NEURON {
	  POINT_PROCESS NaPulse
    USEION na WRITE ina
    RANGE t0, t1, ina
}

UNITS {
    (nA) = (nanoamp)
}

PARAMETER {		:parameters that can be entered when function is called in cell-setup 
	  celsius = 34	    (degC)
    t0 = 1            (ms)           : Start time of pulse
    t1 = 2            (ms)           : End time of pulse
    inabar = -0.0001 (nA)           : Amplitude of natamate pulse
    : inabar = -0.0000 (nA)           : Amplitude of natamate pulse
}

ASSIGNED {                       : parameters needed to solve DE
	  ina   (nA)
}

INITIAL {                        : initialize the following parameter using rates()
	  ina = 0
}

BREAKPOINT {
    if ((t >= t0) && (t <= t1)) {
        ina = inabar
    } else {
        ina = 0
    }
}
