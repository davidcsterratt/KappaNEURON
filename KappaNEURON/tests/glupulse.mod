TITLE Glutamate injection pulse
: Sets Glu current to iglubar between t0 and t1; conductance is 0 otherwise

NEURON {
	  POINT_PROCESS GluPulse
    USEION glu  WRITE iglu VALENCE 1
    NONSPECIFIC_CURRENT i
    RANGE t0, t1, iglu, i
}

UNITS {
    (nA) = (nanoamp)
}

PARAMETER {		:parameters that can be entered when function is called in cell-setup 
	  celsius = 34	    (degC)
    t0 = 1            (ms)           : Start time of pulse
    t1 = 2            (ms)           : End time of pulse
    iglubar = -0.0001 (nA)           : Amplitude of glutamate pulse
    : iglubar = -0.0000 (nA)           : Amplitude of glutamate pulse
}

ASSIGNED {                       : parameters needed to solve DE
    i      (nA)
	  iglu   (nA)
}

INITIAL {                        : initialize the following parameter using rates()
	  iglu = 0
    i = 0
}

BREAKPOINT {
    if ((t >= t0) && (t <= t1)) {
        iglu = iglubar
    } else {
        iglu = 0
    }
    i = -iglu
}
