COMMENT

NMDA syanapse that writes a glutamate flux (which goes to rule-based
simulation) and reads NMDAi (which is product of rule-based
simulation).

ENDCOMMENT

NEURON {
    POINT_PROCESS NmdaSynUrak
    RANGE  e, g, i, b, ica, iGlu
    NONSPECIFIC_CURRENT i
    USEION ca READ cai,cao WRITE ica
    USEION NMDA READ NMDAi,iNMDA VALENCE 0
    USEION Glu  WRITE iGlu VALENCE 0
    GLOBAL total, mg, q10, taurise, taufast, tauslow, taurise_exp, taufast_exp, tauslow_exp, afast, aslow, normfac, T_exp, K0, delta, fracca
}

INCLUDE "units.inc"

PARAMETER {
    : Time constants and afast from Andrasfalvy & Magee 01
    taurise_exp =    6.46 (ms) <1e-9,1e9>    : rise      
    taufast_exp =  252.5  (ms) <1e-9,1e9>    : fast decay
    tauslow_exp = 1660    (ms) <1e-9,1e9>    : slow decay
    afast = 0.61 <0,1>
    e=0	(mV)
    mg	= 1    (mM)		: external magnesium concentration
    fracca= 0.13        : fraction of current that is ca ions; Srupuston &al 95
    z = 2
    celsius = 22	(degC)
    T_exp = 22    (degC)
    q10 = 3       : Hestrin 90
    K0 = 4.1 (mM) : From Spruston &al 95
    delta = 0.8   : From Spruston &al 95
    N_A = 6.02205E23 
}

ASSIGNED {
    v       (mV)
    i       (nA)
    ica	    (nA) 	
    iGlu    (nA)
    iNMDA   (nA)
    g       (uS)
    aslow 
    total   (uS)
    cai     (mM)
    cao     (mM)
    taurise (ms)
    taufast (ms)
    tauslow (ms)
    normfac 
    b
    NMDAi  (mM)
    vol  (micrometer3)
    area (micrometer2)
    diam (micrometer)
    L    (micrometer)
    conv  (/mM)
    fGlu
}

INITIAL {
    L = area/(PI*diam)
    vol = L*PI*(diam/2)^2
    conv = N_A*vol*(1e-18)
}


BREAKPOINT {
    b = mgblock(v)		: b is the block by magnesium at this voltage
    : Single channel conductnace * conversion factor * concentration of NMDAi * block factor
    g = 0.045E-3(uS) *conv * NMDAi * b 
    : It's not possible to create a flux of an uncharged molecule;
    : therefore glutamate has to be created as a current, which feeds
    : into the rule-based simulation.
    iGlu = fGlu*-0.0001
    : The currents need to be corrected for area of spine, assuming
    : that the synapse is on a spine. Also, effect of the glutamate
    : current needs to be cancelled out.
    i =   L/(L + diam/4) * g * (1-fracca) * (v - e) - iGlu
    ica = L/(L + diam/4) * g * fracca     * ghkg(v,cai,cao,z)
    : Should print out a near-integer number of open NMDA gates
    : printf("NMDAi=%g\n", conv * NMDAi)
}

NET_RECEIVE(weight (uS)) {
    if (flag == 0) {
        : If flag==0, an event has been received from a netstim.
        fGlu = 1
        : This sends an event which arrives in 0.1ms with flag = 1.
        net_send(0.1, 1)
    } else {
        : If flag==1, an event has been recived from this mod
        : file. Thus there will have been a 0.1ms long pulse of
        : glutamate.
        fGlu = 0
    }
}

INCLUDE "ghk.inc"

FUNCTION mgblock(v(mV)) {
    TABLE 
    DEPEND mg, K0, delta, celsius
    FROM -140 TO 80 WITH 1000
    : From Spruston &al 95
    mgblock = 1/(1+(mg/K0)*exp(-delta*z*FARADAY*v*(0.001)/R/(celsius+273)))
}



