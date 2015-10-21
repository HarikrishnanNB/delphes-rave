/*
 *  Delphes: a framework for fast simulation of a generic collider experiment
 *  Copyright (C) 2012-2014  Universite catholique de Louvain (UCL), Belgium
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef SecondaryVertexTagging_h
#define SecondaryVertexTagging_h

/** \class SecondaryVertexTagging
 *
 *  Selects candidates from the InputArray according to the efficiency formula.
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "classes/DelphesModule.h"

#include <vector>
#include <utility>
#ifndef __CINT__
#include <unordered_map>
#endif

class TObjArray;
class Candidate;
namespace rave {
  class ConstantMagneticField;
  class VertexFactory;
  class Vertex;
  class FlavorTagFactory;
  class Ellipsoid3D;
  class Track;
}
class RaveConverter;

struct SortedTracks {
  std::vector<std::pair<double, Candidate*> > first;
  std::vector<Candidate*> second;
  std::vector<Candidate*> all; // both vertices, include low pt
};

class SecondaryVertexTagging: public DelphesModule
{
public:

  SecondaryVertexTagging();
  ~SecondaryVertexTagging();

  void Init();
  void Process();
  void Finish();

private:
  Double_t fPtMin;
  Double_t fDeltaR;
  Double_t fIPmax;
  double fBz;			// magnetic field along z
  double fPrimaryVertexPtMin;
  double fPrimaryVertexD0Max;
  double fPrimaryVertexCompatibility;
  double fHLSecVxCompatibility;
  double fMidLevelSecVxCompatibility;

  TIterator *fItTrackInputArray; //!
  TIterator *fItJetInputArray; //!

  const TObjArray *fTrackInputArray; //!
  const TObjArray *fJetInputArray; //!

  TObjArray *fOutputArray; //!

  std::vector<Candidate*> GetTracks(Candidate*);
  // return a pair: first is selected tracks in the jet, second is selected
  // tracks not in the jet
  SortedTracks SelectTracksInJet(
    Candidate*, const std::unordered_map<unsigned, double>& primary_weight);
  rave::Vertex GetPrimaryVertex();
  rave::Vertex getPrimaryVertex(const std::vector<rave::Track>& tracks);

  rave::ConstantMagneticField* fMagneticField;
  rave::VertexFactory* fVertexFactory;
  RaveConverter* fRaveConverter;
  rave::FlavorTagFactory* fFlavorTagFactory;
  rave::Ellipsoid3D* fBeamspot;
  std::map<std::string, int> fDebugCounts;

  ClassDef(SecondaryVertexTagging, 1)
};

#endif // SecondaryVertexTagging_h
