#ifndef FLAVOR_TAG_TRUTH_HH
#define FLAVOR_TAG_TRUTH_HH

#include <ostream>

struct TruthVertex
{
  double x;
  double y;
  double z;
  int pdgid;
  int idx;
  int n_charged_tracks;
};

bool operator<(const TruthVertex&, const TruthVertex&);
std::ostream& operator<<(std::ostream&, const TruthVertex&);
#endif
