COMMENT
	calcium accumulation into a volume of area*depth next to the
	membrane with a decay (time constant tau) to resting level
	given by the global calcium variable cai0_ca_ion
ENDCOMMENT

UNITS {
	(mM) = (milli/liter)
	(mA) = (milliamp)
  PI = (pi) (1)
	F = (faraday) (coulombs)
}

NEURON {
	SUFFIX caBuffer
	USEION ca READ cai, ica WRITE cai
	RANGE cai0, B0, gamma1, gamma2, k1, k2, vol
}

PARAMETER {
	  cai0 = 50e-6 (mM)	: Requires explicit use in INITIAL
		: block for it to take precedence over cai0_ca_ion
		: Do not forget to initialize in hoc if different
		: from this default.
    : B0 = 10000 molecules/(_conversion_factor*volume)
    : = 10000/(602214.129 * 0.0785398163397)
    B0 = 0.2114264 (mM)
    gamma1 = 1 (/ms)
    gamma2 = 1 (/ms)
    N_A = 6.02205E23
}

ASSIGNED {
	  ica (mA/cm2)
    k1 (/mM/ms)
    k2 (/ms)
    vol  (micrometer3)
    area (micrometer2)
    diam (micrometer)
    L    (micrometer)
}

STATE {
    cai (mM)
    B   (mM) 
}

INITIAL {
	  cai = cai0
    B = B0
    vol = L*PI*(diam/2)^2
    k1 = (1e-18)*gamma1*N_A*vol
    k2 = gamma2
}

BREAKPOINT {
	SOLVE integrate METHOD derivimplicit
}

DERIVATIVE integrate {
	  cai' = -ica/diam/PI/F/2 * (1e4) -k1*cai*B
    B' = -k1*cai*B + k2*(B0 - B)
}
