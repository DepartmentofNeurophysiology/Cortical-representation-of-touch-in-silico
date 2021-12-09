NEURON {
:  THREADSAFE
  POINT_PROCESS FluctExp2Syn
  RANGE tau_rise, tau_fall, cn, mean_amp, cv, type, e, i, s, g, std, std0, pf, plas, tau_plas
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
  std0 = 1.0
  pf = 0.0
  plas = 0.75
  tau_plas = 120.0
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
  std
}

INITIAL {
  s = 0
  g = 0
  std = std0
  set_seed(seed)
}

BREAKPOINT {
  SOLVE state METHOD cnexp
  : obtaining current, whether it is excitatory (default), denoted by type = 1 (exc), or inhibitory, denoted by type = -1
  if (type == -1) {  : inhibitory
    i = -g               : current positive (g is negative according to the incoming weights in NET_RECIVE
  } else {
    i = g * (v-e) / 70    : current negative
    }
}

DERIVATIVE state {
  s' = -s/tau_rise
  g' = (cn*s-g)/tau_fall
  std' = -(std-1.0)/tau_plas
}

NET_RECEIVE(w (uS)) {
  LOCAL ww, pfail, prob

  : w should be equal to mean_amp
  : stochastic amplitude
  ww = fabs(w*(1+cv*normrand(0,1)))
  : clipping response
  if (ww > 3*fabs(w)){
    ww = 3*fabs(w)
  }
  : add sign
  if (w<0) {
    ww = -ww
  }

  : Includes short-term plasticity
  ww = std * ww

  : failure of transmission
  pfail = pf / std
  
  :scop_random uniform between 0 and 1
  prob = scop_random()
  if (prob >= pfail) {

    s = s + ww         : update conductance

    std = std * plas   : update short-term dynamics
    : clipping
    if (std >= 5) {
      std = 5
    }
    if (std <= 0.4) {
      std = 0.4
    }
  } 
}


COMMENT
with printing options (to check everything is going fine)

NET_RECEIVE(w (uS)) {
  LOCAL ww, pfail, prob

  : w should be equal to mean_amp
  if (flag_print == 1){
    VERBATIM
      printf("Mean amplitude: %g - ", mean_amp);
    ENDVERBATIM
  }

  : stochastic amplitude
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
      printf("Instantiated Amplitude: %g -> ", _lww);
    ENDVERBATIM
  }

  : Includes short-term plasticity
  ww = std * ww
  if (flag_print == 1){
    VERBATIM
      printf("Potentiated: %g\n", _lww);
    ENDVERBATIM
  }

  : failure of transmission
  pfail = pf / std
  
  :scop_random uniform between 0 and 1
  prob = scop_random()

  if (flag_print == 1){
    VERBATIM
      printf(" Failure probability: %g (prob: %g)\n", _lpfail, _lprob);
      printf("  Conductance before (putative) event: %g\n", s);
    ENDVERBATIM
  }

  if (prob >= pfail) {

    s = s + ww         : update conductance

    if (flag_print == 1){
      VERBATIM
        printf("    Event transmitted\n");
        printf("       STD Plasticity before: %g\n", std);
      ENDVERBATIM
    }
    std = std * plas   : update short-term dynamics
    : clipping
    if (std >= 5) {
      std = 5
    }
    if (std <= 0.4) {
      std = 0.4
    }
  
    if (flag_print == 1){
      VERBATIM
        printf("       STD Plasticity after: %g\n", std);
      ENDVERBATIM
    }

  } else {
    if (flag_print == 1){
      VERBATIM
        printf("    Event not transmitted\n");
      ENDVERBATIM
    }
  }

  if (flag_print == 1){
    VERBATIM
      printf("  Conductance after (putative) event: %g\n", s);
    ENDVERBATIM
  }
 
}

ENDCOMMENT
