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

#ifndef SecondaryVertexAssociator_h
#define SecondaryVertexAssociator_h

/** \class SecondaryVertexAssociator
 *
 *  Find origin of jet and evaluate jet flavor
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "classes/DelphesModule.h"
#include "classes/DelphesClasses.h"
#include "classes/flavortag/flavor_tag_truth.hh"

#include <map>
#include <set>

class TObjArray;


class SecondaryVertexAssociator: public DelphesModule
{
public:

  SecondaryVertexAssociator();
  ~SecondaryVertexAssociator();

  void Init();
  void Process();
  void Finish();

private:

  TIterator *fItJetInputArray; //!

  const TObjArray *fParticleInputArray; //!
  const TObjArray *fJetInputArray; //!

  typedef std::vector<TruthVertex> HFVs;
  HFVs getHeavyFlavorVertices(Candidate* track);
  HFVs walkTruthRecord(Candidate* genPart, const std::set<int>& targets);
  std::vector<Candidate*> getStableChildren(Candidate* idx);
  Candidate* getGenPart(int idx);

  ClassDef(SecondaryVertexAssociator, 1)
};



#endif
