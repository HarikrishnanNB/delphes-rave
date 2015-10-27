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


/** \class TrackBasedBTagging
 *
 *  b-tagging algorithm based on counting tracks with large impact parameter
 *
 *  \author M. Selvaggi - UCL, Louvain-la-Neuve
 *
 */

#include "modules/TrackBasedBTagging.h"
#include "classes/flavortag/hl_vars.hh"

#include "classes/DelphesClasses.h"
#include "classes/DelphesFactory.h"
#include "classes/DelphesFormula.h"

#include "TMath.h"
#include "TString.h"
#include "TFormula.h"
#include "TRandom3.h"
#include "TObjArray.h"
#include "TLorentzVector.h"

#include <algorithm>
#include <stdexcept>
// #include <iostream>
#include <sstream>

using namespace std;

//------------------------------------------------------------------------------

TrackBasedBTagging::TrackBasedBTagging() :
  fItTrackInputArray(0), fItJetInputArray(0)
{
}

//------------------------------------------------------------------------------

TrackBasedBTagging::~TrackBasedBTagging()
{
}

//------------------------------------------------------------------------------

void TrackBasedBTagging::Init()
{
  fPtMin = GetDouble("TrackPtMin", 1.0);
  fDeltaR = GetDouble("DeltaR", 0.3);
  fIPmax = GetDouble("TrackIPMax", 2.0);

  // import input array(s)

  fTrackInputArray = ImportArray(GetString("TrackInputArray", "Calorimeter/eflowTracks"));
  fItTrackInputArray = fTrackInputArray->MakeIterator();

  fJetInputArray = ImportArray(GetString("JetInputArray", "FastJetFinder/jets"));
  fItJetInputArray = fJetInputArray->MakeIterator();
}

//------------------------------------------------------------------------------

void TrackBasedBTagging::Finish()
{
  if(fItTrackInputArray) delete fItTrackInputArray;
  if(fItJetInputArray) delete fItJetInputArray;
}

//------------------------------------------------------------------------------

void TrackBasedBTagging::Process()
{

  // loop over all input jets
  fItJetInputArray->Reset();
  Candidate* jet;
  while((jet = static_cast<Candidate*>(fItJetInputArray->Next())))
  {
    const TLorentzVector &jetMomentum = jet->Momentum;

    std::vector<TrackParameters> trk_pars;
    if (jet->GetTracks()->GetEntriesFast() > 0) {
      throw std::logic_error("tried to add traks to a jet twice");
    }
    // loop over all input tracks
    fItTrackInputArray->Reset();
    Candidate* track;
    while((track = static_cast<Candidate*>(fItTrackInputArray->Next())))
    {
      const TLorentzVector &trkMomentum = track->Momentum;

      double dr = jetMomentum.DeltaR(trkMomentum);
      double tpt = trkMomentum.Pt();
      double dxy = std::abs(track->Dxy);

      if(tpt < fPtMin) continue;
      if(dr > fDeltaR) continue;
      if(dxy > fIPmax) continue;
      trk_pars.emplace_back(track->trkPar, track->trkCov);
      // std::cout << trk_pars.back() << std::endl;
      jet->AddTrack(track);
    }
    jet->hlTrk.fill(jetMomentum.Vect(), trk_pars);
    // std::cout << jet->hlTrk << std::endl;

  }
}

//------------------------------------------------------------------------------
