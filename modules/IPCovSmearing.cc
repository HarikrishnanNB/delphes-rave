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

/** \class IPCovSmearing
 *
 *  Performs track smearing
 *
 *  \author Shih-Chieh Hsu
 *  \author Dan Guest (bug fixes and Eigen rewrite)
 *
 */


#include "modules/IPCovSmearing.h"

// for TRKCOV_TOARRAY macros
#include "classes/flavortag/track_set_macros.hh"

#include "classes/DelphesClasses.h"
#include "classes/DelphesFactory.h"

#include "ExRootAnalysis/ExRootResult.h"
#include "ExRootAnalysis/ExRootFilter.h"
#include "ExRootAnalysis/ExRootClassifier.h"

#include "TFile.h"
#include "TFormula.h"
#include "TObjArray.h"
#include "TLorentzVector.h"
#include "TMatrixDSym.h"

// Eigen
#include <Eigen/Cholesky>

#include <algorithm>
#include <stdexcept>
#include <iostream>
#include <memory>

namespace {
  const double pi = std::atan2(0, -1);

  // CovMatrix is a typedef for an Eigen matrix, defined in header
  CovMatrix get_cov_matrix(TFile&, int ptbin, int etabin);

  // The smearing matrices don't include a bin below 10 GeV, we hack
  // this instead.
  void do_low_pt_hack(CovMatrix& matrix);

  // Covairance matrices were defined in MeV, need to convert to GeV
  void convert_units_to_gev(CovMatrix&);

  // Function to copy to the Candidate
  void set_covariance(float*, const CovMatrix& matrix);
}

//------------------------------------------------------------------------------

IPCovSmearing::IPCovSmearing() :
  fNBinMisses(0),
  fItInputArray(0)
{
}

//------------------------------------------------------------------------------

IPCovSmearing::~IPCovSmearing()
{
}

//------------------------------------------------------------------------------

void IPCovSmearing::Init()
{

  double smear_mult = GetDouble("SmearingMultiple", 1.0);

  const char* filename = GetString("SmearParamFile", "Parametrisation/IDParametrisierung.root");

  TFile file_para(filename,"READ");
  if (!file_para.IsOpen() || file_para.IsZombie()) {
    throw std::runtime_error("bad file: " + std::string(filename));
  }


  ptbins.push_back(10);
  ptbins.push_back(20);
  ptbins.push_back(50);
  ptbins.push_back(100);
  ptbins.push_back(200);
  ptbins.push_back(250);
  ptbins.push_back(500);
  ptbins.push_back(750);

  etabins.push_back(0.0);
  etabins.push_back(0.4);
  etabins.push_back(0.8);
  etabins.push_back(1.05);
  etabins.push_back(1.5);
  etabins.push_back(1.7);
  etabins.push_back(2.0);
  etabins.push_back(2.25);
  etabins.push_back(2.7);

  const int pt_bin_max = ptbins.size();
  const int eta_bins_max = etabins.size();
  for (int ipt = -1 ; ipt < pt_bin_max; ipt++) {
    for (int ieta = 0; ieta < eta_bins_max; ieta++) {
      try {
	CovMatrix cov = get_cov_matrix(file_para, ipt, ieta) * smear_mult;
	fCovarianceMatrices[ipt][ieta] = cov;

	// get the lower part of the Cholesky decomposition. The smearing
	// will be s = L*r, where r is a random gaussian 5-vector
	fSmearingMatrices[ipt][ieta] = cov.llt().matrixL();

      } catch (std::invalid_argument&) {
	std::cout << "** INFO: no smearing defined for pt-eta "
		  << ipt << " " << ieta << std::endl;
      }
    }
  }

  // import input array

  fInputArray = ImportArray(GetString("InputArray", "TrackMerger/tracks"));
  fItInputArray = fInputArray->MakeIterator();

  // create output array

  fOutputArray = ExportArray(GetString("OutputArray", "tracks"));
}


//------------------------------------------------------------------------------

void IPCovSmearing::Finish()
{
  if(fItInputArray) delete fItInputArray;
  if (fNBinMisses) {
    std::cout << "PROBLEM: " << fNBinMisses << "bin misses in track smearing"
	      << std::endl;
  }
}

//------------------------------------------------------------------------------


void IPCovSmearing::Process()
{
  Candidate *track;
  double xd, yd, zd;
  double pt, eta, px, py, phi;
  int charge;

  fItInputArray->Reset();
  while((track = static_cast<Candidate*>(fItInputArray->Next())))
  {
    // Delphes gives us a ``track'' candidate, which has a single
    // sub-candidate, the original generated particle. We smear the
    // track, and insert the original track _between_ the new track
    // and the generated particle.
    //
    // We take momentum and position from the generated particle to
    // avoid applying smearing twice.
    Candidate* particle = static_cast<Candidate*>(
      track->GetCandidates()->At(0));

    // check the above assumption: the particle should have no
    // sub-candidates.
    assert(particle->GetCandidates()->GetEntriesFast() == 0);

    const TLorentzVector &candidateMomentum = particle->Momentum;

    charge = particle->Charge;

    eta = candidateMomentum.Eta();
    pt = candidateMomentum.Pt();

    // NOTE: the phi used here isn't _strictly_ correct, since it
    //       doesn't extrapolate all the way to perigee.  The actual
    //       measured phi would be deflected slightly more by the
    //       magnetic field. We're using it for consistency with the
    //       rest of Delphes, though.
    phi = candidateMomentum.Phi();

    // Same logic as above should follow here. Again, we'll stick with
    // this for consistency with the rest of delphes (since this
    // should be a small effect).
    px = candidateMomentum.Px();
    py = candidateMomentum.Py();

    // The d0 and z0 parameters aren't stored in the Gen Particle,
    // they have to be taken from the track. These parameters _are_
    // taken from perigee, unlike px and py above.
    xd =  track->Xd;
    yd =  track->Yd;
    zd =  track->Zd;

    double phid0 = phi - pi/2;

    // Compute qoverp and theta: Because matrix parametrisation is for
    // (d0,z0,phi,theta,qoverp)
    double qoverp = charge/(pt*cosh(eta));
    double theta = 2.*std::atan(std::exp(-eta));

    // calculate impact parameter (_before_ smearing). Note that this
    // isn't the true perigee point (see notes above).
    double d0 = (xd*py - yd*px)/pt;
    double z0 = zd;

    // get pt and eta bins (TODO: replace with something less confusing)
    int ptbin = -1;
    for(unsigned int i=0;i< ptbins.size();i++){
      if(pt > ptbins.at(i)) ptbin=i;
    }
    int etabin = -1;
    for(unsigned int i=0;i< etabins.size();i++){
      if(fabs(eta) > etabins.at(i)) etabin=i;
    }

    // Now do the smearing
    const auto& bins = getValidBins(ptbin, etabin);
    const CovMatrix& smearing_matrix =
      fSmearingMatrices.at(bins.first).at(bins.second);
    TrackVector track_parameters;
    track_parameters << d0, z0, phi, theta, qoverp;
    TrackVector rand = getRandomVector();
    // calculate the smeared track
    TrackVector smearing = smearing_matrix * rand;
    TrackVector smeared = smearing + track_parameters;

    // Save the current particle as the unsmeared one, then clone and
    // copy smeared parameters to the particle.
    Candidate* smeared_track = static_cast<Candidate*>(track->Clone());

    // copy track parameters to the track
    float* trkPar = smeared_track->trkPar;
    for (int iii = 0; iii < 5; iii++) trkPar[iii] = smeared(iii);
    float* cov_array = smeared_track->trkCov;

    // copy covariance matrix to the track
    const CovMatrix& cov = fCovarianceMatrices.at(bins.first).at(bins.second);
    set_covariance(cov_array, cov);

    // fill the track parameters
    {
      using namespace TrackParam;
      // assign track parameters to the track momentum
      double smeared_pt = charge/(smeared(QOVERP)*cosh(eta));
      assert(smeared_pt >= 0);
      double smeared_eta = -std::log(std::tan(smeared(THETA)/2));
      smeared_track->Momentum.SetPtEtaPhiM(
      smeared_pt, smeared_eta, smeared(PHI), candidateMomentum.M());

      // assign the IP smearing
      double smeared_d0 = smeared(D0);
      smeared_track->Dxy = smeared_d0;
      smeared_track->SDxy = std::sqrt(std::abs(cov_array[D0D0]));

      // smear the Xd and Yd consistent with d0 smearing
      double phid0_reco = phid0 + (smeared(PHI) - phi);
      smeared_track->Xd = smeared_d0 * std::cos(phid0_reco);
      smeared_track->Yd = smeared_d0 * std::sin(phid0_reco);
      smeared_track->Zd = smeared(Z0);
    }

    // Remove the previous gen particle from candidates, insert the
    // original track in the candidate array. The resulting structure
    // is: Smeared Track -> Unsmeared Track -> Gen Particle.
    auto* array = static_cast<TObjArray*>(smeared_track->GetCandidates());
    array->Clear() ;
    array->Add(track);

    fOutputArray->Add(smeared_track);
  }
}

std::pair<int,int> IPCovSmearing::getValidBins(int ptbin, int etabin) {
  const auto& pt_mats = fSmearingMatrices[ptbin];
  // use the lower pt bin if this one isn't defined
  if (!pt_mats.count(etabin)) {
    if (etabin == 0) {
      throw std::logic_error(
	"no eta bins for pt bin: " + std::to_string(etabin));
    } else {
      fNBinMisses++;
      return getValidBins(ptbin, etabin-1);
    }
  }
  return {ptbin, etabin};
}

TrackVector IPCovSmearing::getRandomVector() {
  std::normal_distribution<double> norm(0,1);
  TrackVector pars;
  for (int iii = 0; iii < 5; iii++) {
    pars(iii) = norm(fRandomGenerator);
  }
  return pars;
}

namespace {
  CovMatrix get_cov_matrix(TFile& file, int ptbin, int etabin) {
    bool lowpt_hack = false;
    if (ptbin == -1) {
      lowpt_hack = true;
      ptbin = 0;
    }
    std::unique_ptr<TMatrixDSym> cov(new TMatrixDSym(5));
    TMatrixDSym* cptr = cov.get();
    TString name;
    name.Form("covmat_ptbin%.2i_etabin%.2i",ptbin,etabin);
    file.GetObject(name,cptr);
    if(!cptr){
      throw std::invalid_argument("no bin " + std::string(name.Data()));
    }
    CovMatrix covariance;
    for (int iii = 0; iii < 5; iii++) {
      for (int jjj = 0; jjj < 5; jjj++) {
	covariance(iii,jjj) = (*cptr)(iii,jjj);
      }
    }

    // various conversions
    convert_units_to_gev(covariance);
    if (lowpt_hack) do_low_pt_hack(covariance);

    return covariance;
  }
  void do_low_pt_hack(CovMatrix& cov_matrix){
    // hack to give larger uncertainty to low pt bins
    CovMatrix hack_matrix = CovMatrix::Identity();
    const double unct_mul = 2.0; // uncertainty increase for low pt

    // add non-1 entries
    using namespace TrackParam;
    for (auto comp: {D0, Z0} ) hack_matrix(comp,comp) = unct_mul;

    cov_matrix = hack_matrix * cov_matrix * hack_matrix;

  }
  void convert_units_to_gev(CovMatrix& matrix) {
    CovMatrix gev_from_mev = CovMatrix::Identity();
    gev_from_mev(4,4) = 1000;
    matrix = gev_from_mev * matrix * gev_from_mev;
  }

  void set_covariance(float* cov_array, const CovMatrix& cov) {
    using namespace TrackParam;
    TRKCOV_2TOARRAY(D0);

    TRKCOV_2TOARRAY(Z0);
    TRKCOV_TOARRAY (Z0,D0);

    TRKCOV_2TOARRAY(PHI);
    TRKCOV_TOARRAY (PHI,D0);
    TRKCOV_TOARRAY (PHI,Z0);

    TRKCOV_2TOARRAY(THETA);
    TRKCOV_TOARRAY (THETA,D0);
    TRKCOV_TOARRAY (THETA,Z0);
    TRKCOV_TOARRAY (THETA,PHI);

    TRKCOV_2TOARRAY(QOVERP);
    TRKCOV_TOARRAY (QOVERP,D0);
    TRKCOV_TOARRAY (QOVERP,Z0);
    TRKCOV_TOARRAY (QOVERP,PHI);
    TRKCOV_TOARRAY (QOVERP,THETA);
  }
}

//------------------------------------------------------------------------------
