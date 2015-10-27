#ifndef CONSTANTS_JETPROB_HH
#define CONSTANTS_JETPROB_HH

// JetProb parameters taken from ATLAS-CONF-2010-091
// to fit the distribution:
// p0 * exp(-x*x / (2 * p1*p1)) + p2 * exp(-x*x / (2*p3*p3))
//    + exp(-p4-p5*abs(x)) + exp(-p6-p7*abs(x))

namespace jetprob {
  const double P0 = 0.32;
  const double P1 = 0.56;
  const double P2 = 0.19;
  const double P3 = 0.89;
  const double P4 = 3.09;
  const double P5 = 0.81;
  const double P6 = 6.07;
  const double P7 = 0.23;
}
#endif
