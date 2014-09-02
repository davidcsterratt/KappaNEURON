COMMENT
  Calcium accumlation into the full volume of the compartment and a
  linear membrane pump.
ENDCOMMENT

UNITS {
	(mM) = (milli/liter)
	(mA) = (milliamp)
	F = (faraday) (coulombs)
}

NEURON {
	SUFFIX caPump1
	USEION ca READ cai, ica WRITE cai, ica
	RANGE cai0, k1
}

PARAMETER {
	  cai0 = 50e-6 (mM)	: Requires explicit use in INITIAL
		: block for it to take precedence over cai0_ca_ion
		: Do not forget to initialize in hoc if different
		: from this default.
    k1 = 47.3 (/ms)
}

ASSIGNED {
	  ica   (mA/cm2)
    diam  (micron)
    ipump (mA/cm2)
}

STATE {
    cai (mM)
}

INITIAL {
	  cai = cai0
}

BREAKPOINT {
	  SOLVE integrate METHOD derivimplicit
    : Add the pump current to ica after the SOLVE
    ica = k1*cai/2*diam*F/(1e4)
}

DERIVATIVE integrate {
    cai' = -2*ica/diam/F*(1e4)
    : We don't write this:
    :   cai' = -2*ica/diam/F*(1e4) - k1*cai
    : because the pump current has been added to ica in the BREAKPOINT
    : block

}
