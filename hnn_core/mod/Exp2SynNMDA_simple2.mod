: Exp2SynNMDA_simple2.mod

NEURON {
    POINT_PROCESS Exp2SynNMDA_simple
    RANGE tau1, tau2, e, i
    NONSPECIFIC_CURRENT i
    RANGE g, mgblock, factor
}

UNITS {
    (nA) = (nanoamp)
    (mV) = (millivolt)
    (uS) = (microsiemens)
}

PARAMETER {
    tau1 = 0.5 (ms) 
    tau2 = 44 (ms) 
    e = 0 (mV)
    mg = 1 (mM)        : external magnesium concentration
}

ASSIGNED {
    v (mV)
    i (nA)
    g (uS)
    factor (1)
}

STATE {
    A (uS)
    B (uS)
}

INITIAL {
    LOCAL tp
    if (tau1/tau2 > 0.9999) {
        tau1 = 0.9999 * tau2
    }
    A = 0
    B = 0
    tp = (tau1 * tau2) / (tau2 - tau1) * log(tau2 / tau1)
    factor = -exp(-tp / tau1) + exp(-tp / tau2)
    factor = 1 / factor
}

BREAKPOINT {
    SOLVE state METHOD cnexp
    g = B - A
    i = g * mgblock(v) * (v - e)
}

DERIVATIVE state {
    A' = -A / tau1
    B' = -B / tau2
}

FUNCTION mgblock(v(mV)) {
    TABLE 
    DEPEND mg
    FROM -140 TO 80 WITH 1000

    : from Jahr & Stevens
    mgblock = 1 / (1 + exp(0.062 (/mV) * -v) * (mg / 3.57 (mM)))
}

NET_RECEIVE(w (uS)) {
    A = A + w * factor
    B = B + w * factor
}
