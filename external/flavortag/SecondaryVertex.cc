#include "SecondaryVertex.hh"

#include "classes/DelphesClasses.h"

SecondaryVertex::SecondaryVertex()
{
  clear();
}
SecondaryVertex::SecondaryVertex(double x, double y, double z):
  TVector3(x, y, z)
{
  clear();
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
}
