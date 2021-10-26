NEURON {
:  THREADSAFE
  POINT_PROCESS FluctExp2Syn
  RANGE tau_rise, tau_fall, cn, mean_amp, cv, type, e, i, s, g
  RANGE seed, flag_print
  NONSPECIFIC_CURRENT i
}

UNITS {
  (nA) = (nanoamp)
  (mV) = (millivolt)
  (uS) = (microsiemens)
}

PARAMETER {
  tau_rise = 1.0 (ms) <1e-9,1e9>
  tau_fall = 2.0 (ms) <1e-9,1e9>
  cn = 4         : (tau_rise/tau_fall)^[-tau_fall/(tau_fall-tau_rise)]
  mean_amp = 0.001
  cv = 0.0
  type = 1
  e=0	(mV)
  seed = 12
  flag_print = 0
}

ASSIGNED {
  v (mV)
  i (nA)
}

STATE {
  s (uS)
  g (uS)
}

INITIAL {
  s = 0
  g = 0
  set_seed(seed)
}

BREAKPOINT {
  SOLVE state METHOD cnexp
  : obtaining current, whether it is excitatory (default), denoted by type = 1 (exc), or inhibitory, denoted by type = -1
  if (type == -1) {  : inhibitory
    i = -g               : current positive (g is negative according to the incoming weights in NET_RECIVE
  }
  else {
   i = g * (v-e) / 70    : current negative
  }
}

DERIVATIVE state {
  s' = -s/tau_rise
  g' = (cn*s-g)/tau_fall
}

NET_RECEIVE(w (uS)) {
  LOCAL ww, ww0

  : w should be equal to mean_amp
  ww0 = w

  ww = fabs(w*(1+cv*normrand(0,1)))
  : clipping response
  if (ww > 3*fabs(w)){
    ww = 3*fabs(w)
  }
  : add sign
  if (w<0) {
    ww = -ww
  }
 
  if (flag_print == 1){
    VERBATIM
      printf("%g - %g - %g\n", mean_amp, _lww0, _lww);
    ENDVERBATIM
  }
 
  s = s + ww
}