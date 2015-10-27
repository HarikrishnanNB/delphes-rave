#include "math.hh"

#include <cmath>
#include <cassert>

namespace {
  const double pi = std::atan2(0, -1);
}

double phi_mpi_pi(double phi1, double phi2) {
  double diff = phi1 - phi2;
  if (std::abs(diff) > pi) {
    diff -= std::copysign(2*pi, diff);
  }
  assert(std::abs(diff) <= pi);
  return diff;
}
