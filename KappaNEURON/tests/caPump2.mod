COMMENT
  Calcium accumlation into the full volume of the compartment and a
	nonlinear membrane pump.
ENDCOMMENT

UNITS {
	(mM) = (milli/liter)
	(mA) = (milliamp)
	F = (faraday) (coulombs)
}

NEURON {
	SUFFIX caPump2
	USEION ca READ cai, ica WRITE cai, ica
	RANGE cai0, P0, k1, k2, ipump, P
}

PARAMETER {
	  cai0 = 50e-6 (mM)	: Requires explicit use in INITIAL
		: block for it to take precedence over cai0_ca_ion
		: Do not forget to initialize in hoc if different
		: from this default.
    : P0 = 10000 molecules/(_conversion_factor*volume)
    : = 10000/(602214.129 * 0.0785398163397)
    P0 = 0.20 (mM)
    k1 = 47.3 (/mM-ms)
    k2 = 1    (/ms)
}

ASSIGNED {
	  ica  (mA/cm2)
    diam  (micron)
    ipump  (mA/cm2)
}

STATE {
    cai (mM)
    P   (mM) 
}

INITIAL {
	  cai = cai0
    P = P0
}

BREAKPOINT {
	  SOLVE integrate METHOD derivimplicit
    ipump = k2*(P0 - P)/2*diam*F/(1e4)
    ica = ipump
}

DERIVATIVE integrate {
    : We have to subtract ipump here, because it has already been
    : added to ica in the BREAKPOINT block.
	  cai' = -2*(ica - ipump)/diam/F*(1e4) - k1*cai*P
    P'   = -k1*cai*P + 2*ipump/diam/F*(1e4)
}
