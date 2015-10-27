#include "flavor_tag_truth.hh"


bool operator<(const TruthVertex& v1, const TruthVertex& v2) {
  return v1.idx < v2.idx;
}


std::ostream& operator<<(std::ostream& out, const TruthVertex& vx) {
  out << "#" << vx.idx;
  out << " (" << vx.x << " " << vx.y << " " << vx.z << ") PID: " << vx.pdgid;
  out << ", " << vx.n_charged_tracks << " charged tracks";
  return out;
}
