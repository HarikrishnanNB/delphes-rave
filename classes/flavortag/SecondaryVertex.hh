#ifndef SECONDARY_VERTEX_HH
#define SECONDARY_VERTEX_HH

#include "TVector3.h"

class Candidate;

struct SecondaryVertexTrack
{
  double weight;
  double d0;
  double z0;
  double d0err;
  double z0err;
  double pt;
  // these are relative to jet axis
  double dphi;
  double deta;
  Candidate* delphes_track;
};

class SecondaryVertex: public TVector3
{
public:
  SecondaryVertex();
  SecondaryVertex(double, double, double);
  double Lxy;
  double Lsig;
  double decayLengthVariance;
  int nTracks;
  double eFrac;
  double mass;
  std::string config;
  double deta;
  double dphi;
  // std::vector<std::pair<double, Candidate*> > tracks;
  std::vector<SecondaryVertexTrack> tracks_along_jet;
  void clear();
};

#endif
