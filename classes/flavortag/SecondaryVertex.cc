#include "SecondaryVertex.hh"

#include "classes/DelphesClasses.h"

SecondaryVertex::SecondaryVertex()
{
  clear();
}
SecondaryVertex::SecondaryVertex(double x, double y, double z)
{
  clear();
  SetXYZ(x, y, z);
}
void SecondaryVertex::clear() {
  SetXYZ(0, 0, 0);
  Lxy = -1;
  Lsig = -1;
  nTracks = -1;
  eFrac = -1;
  mass = -1;
  config = "null";
  tracks.clear();
  tracks_along_jet.clear();
}
