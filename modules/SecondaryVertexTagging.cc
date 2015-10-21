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


/** \class SecondaryVertexTagging
 *
 *  Builds secondary vertices from tracks. Uses Rave.
 *
 *  \author Dan Guest
 *
 */

#include "modules/SecondaryVertexTagging.h"

// only define use what's below if we have Rave
#include "flavortag/math.hh"
#include "flavortag/SecondaryVertex.hh"
#include "classes/DelphesClasses.h"
#include "ExRootAnalysis/ExRootConfReader.h"

#include "TObjArray.h"

#ifndef NO_RAVE 		// check for NO_RAVE flag

#include "rave/Version.h"
#include "rave/Track.h"
#include "rave/Exception.h"
#include "RaveBase/Converters/interface/RaveStreamers.h"
#include "RaveBase/Converters/interface/PerigeeToRaveObjects.h"
#include "RaveBase/Converters/interface/RaveToPerigeeObjects.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "rave/Vector6D.h"
#include "rave/VertexFactory.h"
#include "rave/FlavorTagFactory.h"
#include "rave/ConstantMagneticField.h"
#include "rave/VacuumPropagator.h"

#include <iostream>
#include <map>
#include <set>
#include <iomanip>

// forward declare some utility functions that are used below
namespace {
  // assume the pion hypothisis for all tracks
  const double M_PION = 139.57e-3; // in GeV
  // some vertex properties require that we cut on a track association
  // probibility.
  const double VPROB_THRESHOLD = 0.5;

  typedef std::vector<std::pair<double, Candidate*> > WeightedTracks;

  // - walk up the candidate tree to find the generated particle
  Candidate* get_part(Candidate* cand);
  // - dump info about a track
  void print_track_info(const Candidate* cand);
  void print_more_info(const Candidate* cand);
  void print_rave_track_info(const rave::Track& track);
  // - vertex significance
  double vertex_significance(const rave::Vertex&);
  double decay_length_variance(const rave::Vertex&);

  // Several functions to get discriminating tagging info notes:
  //  - These generally assume the pion hypothisis.
  //  - A rave vertex includes a collection of weighted tracks.
  //    In some cases we multiply the weight by the quantity in quesiton,
  //    while in others we only consider tracks above some threshold.
  //  - TODO: rationalize this a bit. Why use thresholds at all?
  double vertex_energy(const rave::Vertex&);
  double cut_vertex_energy(const rave::Vertex&, double threshold);
  double weighted_vertex_energy(const rave::Vertex&);
  double track_energy(const std::vector<rave::Track>&);
  double track_energy(const std::vector<Candidate*>&);
  int n_tracks(const rave::Vertex&, double threshold);
  double mass(const rave::Vertex&, double threshold);
  WeightedTracks delphes_tracks(const rave::Vertex&);
  std::vector<SecondaryVertexTrack> get_tracks_along_jet(
    const WeightedTracks& tracks, const TVector3& jet, double threshold);
  int get_n_shared(const std::vector<SecondaryVertex>& vertices);

  std::ostream& operator<<(std::ostream& os, const SecondaryVertex&);
  std::ostream& operator<<(std::ostream&, const rave::PerigeeParameters5D&);
  std::string oneline(std::string);
}

// class to convert Delphes stuff to Rave stuff
class RaveConverter
{
public:
  RaveConverter(double Bz, double cov_scaling = 1);
  std::vector<rave::Track> getRaveTracks(const std::vector<Candidate*>& in);
  rave::Point3D getSeed(const std::vector<Candidate*>& in);
private:
  rave::Vector6D getState(const Candidate*);
  rave::PerigeeCovariance5D getPerigeeCov(const Candidate*);
  double getRho(double pt_in_gev, int charge); // return in cm^-1
  double getRhoAlt(double qoverp, double theta);
  double getQOverP(double rho, double theta);
  double getDrhoDqoverp(double theta);
  double getDrhoDtheta(double qoverp, double theta);
  double _bz;
  double _cov_scaling;
  rave::PerigeeToRaveObjects _to_rave;
  rave::RaveToPerigeeObjects _to_perigee;
};


//------------------------------------------------------------------------------

SecondaryVertexTagging::SecondaryVertexTagging() :
  fItTrackInputArray(0), fItJetInputArray(0), fMagneticField(0),
  fVertexFactory(0), fRaveConverter(0), fFlavorTagFactory(0), fBeamspot(0)
{
}

//------------------------------------------------------------------------------

SecondaryVertexTagging::~SecondaryVertexTagging()
{
  delete fMagneticField;
  delete fVertexFactory;
  delete fRaveConverter;
  delete fFlavorTagFactory;
  delete fBeamspot;
}

//------------------------------------------------------------------------------

namespace {
  // you own this pointer, be careful with it
  rave::Ellipsoid3D* new_beamspot(ExRootConfParam beamspot_params,
				  std::vector<double> default_beamspot) {
    std::vector<double> beamspot;
    int npars = beamspot_params.GetSize();
    for (int iii = 0; iii < npars; iii++){
      beamspot.push_back(beamspot_params[iii].GetDouble());
    }
    if (beamspot.size() == 0) {
      beamspot = default_beamspot;
    } else if (beamspot.size() != 3) {
      throw std::runtime_error(
	"Beamspot should be specified by sig_x, sig_y, sig_z");
    }
    using namespace rave;
    using namespace std;
    Point3D point(0,0,0);
    // convert beamspot width to cm, and square for variance
    double xx = pow(beamspot.at(0)*0.1, 2);
    double yy = pow(beamspot.at(1)*0.1, 2);
    double zz = pow(beamspot.at(2)*0.1, 2);
    Covariance3D cov(xx, 0, 0,
		        yy, 0,
		            zz);
    return new Ellipsoid3D(point, cov);
  }
  std::vector<std::string> get_string(ExRootTask* rdr,
				      std::string param){
    std::vector<std::string> outs;
    ExRootConfParam string_params = rdr->GetParam(param.c_str());
    int npars = string_params.GetSize();
    if (npars == 0) {
      throw std::runtime_error("missing param " + param);
    }
    for (int iii = 0; iii < npars; iii++) {
      outs.push_back(string_params[iii].GetString());
    }
    return outs;
  }
  rave::Track get_ghost(const TVector3& jet);
  int n_over(const std::vector<std::pair<float, rave::Track> >& trks,
	     float threshold = 0.5);
  std::string avr_config(double vx_compat);
  std::string avf_config(double vx_compat);
  SecondaryVertex sv_from_rave_sv(const rave::Vertex&, double jet_track_e,
				  const TVector3& jet, double threshold = 0);
}

void SecondaryVertexTagging::Init()
{
  // tracking parameters
  fPtMin = GetDouble("TrackPtMin", 1.0);
  fDeltaR = GetDouble("DeltaR", 0.3);
  fIPmax = GetDouble("TrackIPMax", 2.0);

  // magnetic field
  fBz = GetDouble("Bz", 2.0);

  // beamspot (should be specified in mm, converted to cm internally)
  const double bs_xy = 15e-3; 	// 15 microns
  fBeamspot = new_beamspot(GetParam("Beamspot"), {bs_xy, bs_xy, 46.0});

  // primary vertex definition
  fPrimaryVertexPtMin = GetDouble("PrimaryVertexPtMin", 1);
  fPrimaryVertexD0Max = GetDouble("PrimaryVertexD0Max", 0.1);
  fPrimaryVertexCompatibility = GetDouble("PrimaryVertexCompatibility", 0.5);
  // rave method
  fHLSecVxCompatibility = GetDouble("HLSecVxCompatibility", 3.0);
  fMidLevelSecVxCompatibility = GetDouble("MidLevelSecVxCompatibility", 1.0);

  // import input array(s)
  fTrackInputArray = ImportArray(
    GetString("TrackInputArray", "Calorimeter/eflowTracks"));
  fItTrackInputArray = fTrackInputArray->MakeIterator();
  fJetInputArray = ImportArray(
    GetString("JetInputArray", "FastJetFinder/jets"));
  fItJetInputArray = fJetInputArray->MakeIterator();

  // create output array(s)

  // TODO: make this actually do something...
  fOutputArray = ExportArray(GetString("OutputArray", "secondaryVertices"));

  // initalize Rave
  std::cout << "** INFO:     This is Rave Version " << rave::Version()
	    << std::endl;
  fMagneticField = new rave::ConstantMagneticField(0, 0, fBz);
  fVertexFactory = new rave::VertexFactory(
    *fMagneticField, rave::VacuumPropagator(), *fBeamspot, "default", 0);

  double cov_scaling = GetDouble("CovarianceScaling", 1.0);
  fRaveConverter = new RaveConverter(fBz, cov_scaling);
  // to do list
  std::cout << "** TODO: - make a b-tagger that works in Delphes\n"
	    << "         ? figure out why we sometimes get NaN for Lsig\n"
	    << "         ? compute mahalanobis distance for rave coords\n"
	    << "         ? quantify delphes Dxy vs actual Dxy\n"
	    << "         ? dig into FlavorTagFactory, see if I can use it\n"
	    << std::flush;
  // edm::setLogLevel(edm::Error);
}

//------------------------------------------------------------------------------

void SecondaryVertexTagging::Finish()
{
  if(fItTrackInputArray) delete fItTrackInputArray;
  if(fItJetInputArray) delete fItJetInputArray;
  if (fDebugCounts.size() != 0) {
    std::cout << std::endl;
    std::cout << "################################" << std::endl;
    std::cout << "##### some things not good #####" << std::endl;
    std::cout << "################################" << std::endl;
  }
  for (const auto& prob: fDebugCounts) {
    std::cout << prob.first << ": " << prob.second << std::endl;
  }
  if (fDebugCounts.size() != 0) std::cout << std::endl;
}

//------------------------------------------------------------------------------


std::vector<Candidate*> SecondaryVertexTagging::GetTracks(Candidate* jet) {
  // loop over all input jets
  std::vector<Candidate*> jet_tracks;

  const TLorentzVector &jetMomentum = jet->Momentum;

  // loop over all input tracks
  fItTrackInputArray->Reset();
  Candidate* track;
  while((track = static_cast<Candidate*>(fItTrackInputArray->Next())))
  {
    const TLorentzVector &trkMomentum = track->Momentum;

    double dr = jetMomentum.DeltaR(trkMomentum);

    double tpt = trkMomentum.Pt();
    double dxy = std::abs(track->Dxy);
    // double ddxy = track->SDxy;

    if(tpt < fPtMin) continue;
    if(dr > fDeltaR) continue;
    if(dxy > fIPmax) continue;

    jet_tracks.push_back(track);
  }
  return jet_tracks;
}

SortedTracks SecondaryVertexTagging::SelectTracksInJet(
  Candidate* jet, const std::unordered_map<unsigned, double>& primary_wts) {
  // loop over all input jets
  SortedTracks tracks;

  const TLorentzVector &jetMomentum = jet->Momentum;

  // loop over all input tracks
  fItTrackInputArray->Reset();
  Candidate* track;
  while((track = static_cast<Candidate*>(fItTrackInputArray->Next())))
  {
    const TLorentzVector &trkMomentum = track->Momentum;

    double dr = jetMomentum.DeltaR(trkMomentum);

    double tpt = trkMomentum.Pt();
    double dxy = std::abs(track->Dxy);
    // double ddxy = track->SDxy;

    if(dxy > fIPmax) continue;
    bool over_pt_threshold = (tpt >= fPtMin);
    bool track_in_jet = (dr <= fDeltaR);

    unsigned tid = track->GetUniqueID();
    double primary_wt = primary_wts.count(tid) ? primary_wts.at(tid) : -1.0;
    bool track_in_primary = (primary_wt > fPrimaryVertexCompatibility);

    if (track_in_jet) {
      tracks.all.push_back(track);
      if (over_pt_threshold) {
	if(track_in_primary) {
	  tracks.first.emplace_back(primary_wt, track);
	} else {
	  tracks.second.push_back(track);
	}
      }
    }
  }
  assert(tracks.all.size() >= (tracks.first.size() + tracks.second.size()));
  return tracks;
}

void SecondaryVertexTagging::Process()
{
  fItJetInputArray->Reset();
  Candidate* jet;
  const auto& primary = GetPrimaryVertex();
  const auto& primary_tracks = primary.weightedTracks();
  if (n_over(primary_tracks, fPrimaryVertexCompatibility) == 0) {
    fDebugCounts["no primary tracks over threshold"]++;
  }
  std::unordered_map<unsigned, double> primary_weight;
  for (const auto& prim: primary_tracks) {
    const auto* cand = static_cast<Candidate*>(
      prim.second.originalObject());
    primary_weight.emplace(cand->GetUniqueID(), prim.first);
  }

  while((jet = static_cast<Candidate*>(fItJetInputArray->Next())))
  {
    const TLorentzVector& jvec = jet->Momentum;

    auto all_tracks = SelectTracksInJet(jet, primary_weight);
    jet->primaryVertexTracks = get_tracks_along_jet(
      all_tracks.first, jvec.Vect(), fPrimaryVertexCompatibility);
    auto jet_tracks = fRaveConverter->getRaveTracks(all_tracks.second);
    double jet_track_energy = track_energy(all_tracks.all);
    assert(jet_track_energy >= track_energy(all_tracks.second));
    // try out methods:
    // - "kalman" only ever makes one vertex
    // - "mvf" crashes...
    // - "avf": Gives each track a `weight', but forms only one vertex.
    // - "avr": Same as AVF except that it forms new vertices when track
    //          weight drops below 50%.
    // - "tkvf": trimmed Kalman fitter
    // rave::Point3D seed = fRaveConverter->getSeed(all_tracks.first);
    // doesn't make sense to get vertices with < 2 tracks...
    std::vector<SecondaryVertex> hl_svx;
    if (jet_tracks.size() >= 2) {
      auto hl_config = avf_config(fHLSecVxCompatibility);
      auto ml_config = avr_config(fMidLevelSecVxCompatibility);
      try {

	// start with high level variables
	auto hl_vert = fVertexFactory->create(jet_tracks, hl_config);
	for (const auto& vert: hl_vert) {
	  auto out_vert = sv_from_rave_sv(
	    vert, jet_track_energy, jvec.Vect(), VPROB_THRESHOLD);
	  out_vert.config = "high-level";
	  hl_svx.push_back(out_vert);
	} // end vertex filling

	// now fill med level
	auto ml_vert = fVertexFactory->create(jet_tracks, ml_config);
	for (const auto& vert: ml_vert) {
	  auto out_vert = sv_from_rave_sv(
	    vert, jet_track_energy, jvec.Vect());
	  out_vert.config = "med-level";
	  jet->secondaryVertices.push_back(out_vert);
	}
      } catch (cms::Exception& e) {
	fDebugCounts[oneline(e.what())]++;
      }
    }	  // end check for two tracks
    // high level (one fitted vertex)
    assert(hl_svx.size() <= 1);
    jet->hlSecVxTracks.clear();
    if (hl_svx.size() > 0) {
      jet->hlSecVxTracks = hl_svx.at(0).tracks_along_jet;
    }
    jet->hlSvx.fill(jvec.Vect(), hl_svx, 0);

    // medium level (multiple vertices)
    jet->mlSvx.fill(jvec.Vect(), jet->secondaryVertices, 0);
  }   // end jet loop
}

rave::Vertex SecondaryVertexTagging::GetPrimaryVertex() {
  // loop over all input tracks
  fItTrackInputArray->Reset();
  Candidate* track;
  std::vector<Candidate*> vxp_tracks;
  while((track = static_cast<Candidate*>(fItTrackInputArray->Next())))
  {
    const TLorentzVector &trkMomentum = track->Momentum;

    if (trkMomentum.Pt() < fPrimaryVertexPtMin) continue;
    if (std::abs(track->Dxy) > fPrimaryVertexD0Max) continue;
    vxp_tracks.push_back(track);
  }

  auto rave_tracks = fRaveConverter->getRaveTracks(vxp_tracks);

  try {
    return getPrimaryVertex(rave_tracks);
  } catch (cms::Exception& e) {
    fDebugCounts[oneline(e.what())]++;
  }
  return rave::Vertex();
}
rave::Vertex SecondaryVertexTagging::getPrimaryVertex(
  const std::vector<rave::Track>& rave_tracks)
{
  auto vertices = fVertexFactory->create(rave_tracks, "avf", true);
  if (vertices.size() == 0) {
    fDebugCounts["no primary vertex"]++;
    return rave::Vertex();
  } else if (vertices.size() > 1) {
    printf("found %lu vertices!?\n", vertices.size());
  }
  const auto& vx = vertices.at(0);
  return vx;
}

namespace {
  std::string avr_config(double vx_compat) {
    std::string vxc = std::to_string(vx_compat);
    return "avr-primcut:" + vxc + "-seccut:" + vxc;
  }
  std::string avf_config(double vx_compat) {
    std::string vxc = std::to_string(vx_compat);
    return "avf-sigmacut:" + vxc;
  }
  SecondaryVertex sv_from_rave_sv(const rave::Vertex& vert,
				  double jet_track_energy,
				  const TVector3& jet,
				  double threshold) {
    auto pos_mm = vert.position() * 10; // convert to mm
    SecondaryVertex out_vert(pos_mm.x(), pos_mm.y(), pos_mm.z());
    out_vert.Lsig = vertex_significance(vert);
    out_vert.Lxy = pos_mm.perp();
    out_vert.decayLengthVariance = decay_length_variance(vert);
    out_vert.nTracks = n_tracks(vert, threshold);
    out_vert.eFrac = cut_vertex_energy(vert, threshold) /
      jet_track_energy;
    out_vert.mass = mass(vert, threshold);
    double vertex_phi = std::atan2(pos_mm.y(), pos_mm.x());
    out_vert.dphi = phi_mpi_pi(vertex_phi, jet.Phi());
    out_vert.deta = out_vert.Eta() - jet.Eta();

    out_vert.tracks = delphes_tracks(vert);
    out_vert.tracks_along_jet = get_tracks_along_jet(
      out_vert.tracks, jet, threshold);
    return out_vert;
  }
  int get_n_shared(const std::vector<SecondaryVertex>& vertices) {
    std::set<Candidate*> tracks;
    std::set<Candidate*> shared;
    for (const auto& vx: vertices) {
      for (const auto& wt_trk: vx.tracks) {
	if (tracks.count(wt_trk.second) ) {
	  shared.insert(wt_trk.second);
	}
	tracks.insert(wt_trk.second);
      }
    }
    return shared.size();
  }

  int n_over(const std::vector<std::pair<float, rave::Track> >& trks,
	     float threshold)
  {
    int n_trk = 0;
    for (const auto& tk: trks) {
      if (tk.first > threshold) n_trk++;
    }
    return n_trk;
  }
  rave::Track get_ghost(const TVector3& jvec) {
    rave::Vector3D rave_jet_momentum(jvec.Px(), jvec.Py(), jvec.Pz());
    rave::Vector6D rave_jet(rave::Point3D(0,0,0), rave_jet_momentum);
    rave::Covariance6D dummy_cov(
      0,0,0,  0,0,   0,	    // position (xx, xy, xz, yy, yz, zz)
      0,0,0,  0,0,0, 0,0,0, // x, y, z vs px, py, pz (dxpx, dxpy, ...)
      0,0,0,  0,0,   0);    // momentum (pxpx, pxpy, ....)
    return rave::Track(rave_jet, dummy_cov, 0, 0, 0);
  }
  // various functions to work with rave (forward declared above)
  double vertex_significance(const rave::Vertex& vx) {
    double Lx = vx.position().x();
    double Ly = vx.position().y();
    double Lz = vx.position().z();
    double decaylength = std::sqrt(Lx*Lx + Ly*Ly + Lz*Lz);
    if (decaylength == 0) return 0;

    double err = sqrt(decay_length_variance(vx));
    return decaylength / err;
  }
  double decay_length_variance(const rave::Vertex& vx) {
    double Lx = vx.position().x();
    double Ly = vx.position().y();
    double Lz = vx.position().z();
    double decaylength = std::sqrt(Lx*Lx + Ly*Ly + Lz*Lz);
    if (decaylength == 0) return 0;

    const rave::Covariance3D& cov = vx.error();
    double xhat = Lx/decaylength;
    double yhat = Ly/decaylength;
    double zhat = Lz/decaylength;
    double var =
      xhat*xhat*cov.dxx() +
      yhat*yhat*cov.dyy() +
      zhat*zhat*cov.dzz() +
      2.*xhat*yhat*cov.dxy() +
      2.*xhat*zhat*cov.dxz() +
      2.*yhat*zhat*cov.dyz();
    return var;
  }

  double vertex_energy(const rave::Vertex& vx) {
    using namespace std;
    double energy = 0;
    for (const auto& trk: vx.tracks()) {
      double tk_e = sqrt(trk.momentum().mag2() + pow(M_PION, 2));
      energy += tk_e;
    }
    return energy;
  }
  double cut_vertex_energy(const rave::Vertex& vx, double threshold) {
    using namespace std;
    double energy = 0;
    for (const auto& wt_trk: vx.weightedTracks()) {
      if (wt_trk.first > threshold) {
	energy += sqrt(wt_trk.second.momentum().mag2() + pow(M_PION, 2));
      }
    }
    return energy;
  }
  double weighted_vertex_energy(const rave::Vertex& vx) {
    // return the weightd track energy: multiply the energy by the
    // track weight in the vertex fit
    using namespace std;
    double energy = 0;
    for (const auto& wt_trk: vx.weightedTracks()) {
      double tk_e = sqrt(wt_trk.second.momentum().mag2() + pow(M_PION, 2));
      energy += tk_e*wt_trk.first;
    }
    return energy;
  }
  double track_energy(const std::vector<rave::Track>& tracks) {
    using namespace std;
    double energy = 0;
    for (const auto& tk: tracks) {
      energy += sqrt(tk.momentum().mag2() + pow(M_PION, 2));
    }
    return energy;
  }
  double track_energy(const std::vector<Candidate*>& tracks) {
    using namespace std;
    double energy = 0;
    for (const auto& tk: tracks) {
      energy += sqrt(tk->Momentum.Vect().Mag2() + pow(M_PION, 2));
    }
    return energy;
  }
  int n_tracks(const rave::Vertex& vx, double threshold) {
    int n_tracks = 0;
    for (const auto& wt_trk: vx.weightedTracks()) {
      if (wt_trk.first > threshold) n_tracks++;
    }
    return n_tracks;
  }
  double mass(const rave::Vertex& vx, double threshold) {
    using namespace std;
    rave::Point3D sum_momentum(0,0,0);
    double sum_energy = 0;
    for (const auto& wt_trk: vx.weightedTracks()) {
      if (wt_trk.first > threshold) {
	const rave::Vector3D& mom = wt_trk.second.momentum();
	sum_momentum += mom;
	sum_energy += sqrt(mom.mag2() + pow(M_PION,2));
      }
    }
    return sqrt(pow(sum_energy,2) - sum_momentum.mag2() );
  }
  WeightedTracks delphes_tracks(const rave::Vertex& vx) {
    WeightedTracks out;
    for (const auto& tkpair: vx.weightedTracks()) {
      auto* cand = static_cast<Candidate*>(tkpair.second.originalObject());
      out.emplace_back(tkpair.first, cand);
    }
    return out;
  }
  std::vector<SecondaryVertexTrack> get_tracks_along_jet(
    const WeightedTracks& delphes_tracks,
    const TVector3& jet, double threshold){
    std::vector<SecondaryVertexTrack> sv_trk;
    for (const auto& wt_trk: delphes_tracks) {
      if (wt_trk.first < threshold) continue;
      const auto& trk = wt_trk.second;
      TrackParameters params(trk->trkPar, trk->trkCov);
      SecondaryVertexTrack track;
      track.weight = wt_trk.first;
      track.d0 = params.d0;
      track.z0 = params.z0;
      track.d0err = params.d0err;
      track.z0err = params.z0err;
      track.pt = trk->Momentum.Pt();
      track.dphi = phi_mpi_pi(params.phi, jet.Phi());
      track.deta = trk->Momentum.Eta() - jet.Eta();
      sv_trk.push_back(track);
    }
    return sv_trk;
  }

  std::ostream& operator<<(std::ostream& os, const SecondaryVertex& vx){
    using namespace std;
    os << " (" << vx.X() << " " << vx.Y() << " " << vx.Z() << ") ";
    os << "eta: " << vx.Eta() << " ";
    os << fixed << right << setprecision(4);
    os << " lxy: " << vx.Lxy;
    os << " sig3d: " << setprecision(1) << setw(5) << vx.Lsig;
    os << " ntrack: " << setw(2) << vx.nTracks;
    os << " efrac: " << setprecision(4) << vx.eFrac;
    os << " mass: " << setw(2) << vx.mass;
    return os;
  }
  std::string oneline(std::string prob){
    std::replace(prob.begin(), prob.end(), '\n','%');
    return prob;
  }
}


//------------------------------------------------------------------------------
// definition of RaveConverter

namespace {
// constants copied from ParticlePropagator
  const double c_light = 2.99792458E8; // in [m/sec]

}

RaveConverter::RaveConverter(double Bz, double cov_scaling):
  _bz(Bz), _cov_scaling(cov_scaling)
{
}

rave::Vector6D RaveConverter::getState(const Candidate* cand) {
  using namespace TrackParam;
  double a_d0 = cand->trkPar[D0];
  double a_z0 = cand->trkPar[Z0];
  double a_phi = cand->trkPar[PHI];
  double a_qoverp = cand->trkPar[QOVERP];
  double a_theta = cand->trkPar[THETA];
  double a_q = cand->Charge;

  // -- translate these to Rave coordinates
  // rave base units are cm and GeV, Delphes takes mm and GeV
  double r_rho = getRhoAlt(a_qoverp, a_theta);
  double r_theta = a_theta;
  double r_phip = a_phi;

  // TODO: check the Delphes Dxy / D0 definition.
  // d0 is x cross p, where x is at perigee, and p is the initial
  // particle momentum. This isn't strictly accurate but it _should_
  // be a small effect for high pT tracks.

  // convert d0 and z0 to cm
  double r_epsilon = a_d0 * 0.1;
  double r_zp = a_z0 * 0.1;

  // build the parameters
  rave::PerigeeParameters5D pars(r_rho, r_theta, r_phip, r_epsilon, r_zp);
  rave::Point3D referencePoint(0, 0, 0); // what is this?
  return _to_rave.convert(pars, a_q, referencePoint);
}

rave::Point3D RaveConverter::getSeed(const std::vector<Candidate*>& trks) {
  rave::Point3D max(0,0,0);
  for (auto& trk: trks) {
    auto* raw = static_cast<Candidate*>(trk->GetCandidates()->At(0));
    auto* particle = static_cast<Candidate*>(raw->GetCandidates()->At(0));
    const auto pos = particle->Position * 0.1;
    rave::Point3D origin(pos.X(), pos.Y(), pos.Z());
    if (origin.mag() > max.mag()) max = origin;
  }
  return max;
}

rave::PerigeeCovariance5D RaveConverter::getPerigeeCov(const Candidate* cand) {
  using namespace TrackParam;
  // -- translate to Rave coordinates
  // rave base units are cm and GeV, Delphes takes mm and GeV
  // need to calculate some things for the jacobian
  const float* par = cand->trkPar;
  double drdq = getDrhoDqoverp(par[THETA]);
  double drdt = getDrhoDtheta(par[QOVERP], par[THETA]);

  const float* cov0 = cand->trkCov;
  float cov[15];
  for (size_t iii = 0; iii < 15; iii++) {
    cov[iii] = cov0[iii] * _cov_scaling;
  }
  float qq = cov[QOVERPQOVERP];
  float qt = cov[QOVERPTHETA];
  float tt = cov[THETATHETA];

  // now multiply out the J*cov*J^T
  float drr = drdq*drdq*qq + drdt*drdt*tt + 2*drdt*drdq*qt;
  float drt = drdt*tt + drdq*qt;

  float drp = cov[QOVERPPHI]*drdq + cov[THETAPHI]*drdt;
  float dtt = cov[THETATHETA];
  float dtp = cov[THETAPHI];
  float dpp = cov[PHIPHI];
  rave::PerigeeCovariance3D cov3d(drr, drt, drp, dtt, dtp, dpp);

  // now the remaining terms. lengths need to be converted to cm
  float dre = (cov[QOVERPD0]*drdq + cov[THETAD0]*drdt) * 0.1;
  float drz = (cov[QOVERPZ0]*drdq + cov[THETAZ0]*drdt) * 0.1;
  float dte = cov[THETAD0]  * 0.1;
  float dtz = cov[THETAZ0]  * 0.1;
  float dpe = cov[PHID0]    * 0.1;
  float dpz = cov[PHIZ0]    * 0.1;
  float dee = cov[D0D0]     * 0.01;
  float dez = cov[Z0D0]     * 0.01;
  float dzz = cov[Z0Z0]     * 0.01;
  rave::PerigeeCovariance5D cov5d(cov3d,
				  dre, drz,
				  dte, dtz,
				  dpe, dpz,
				  dee, dez, dzz);
  return cov5d;
}


double RaveConverter::getRho(double pt_in_gev, int charge) {
  // rave wants rho = 1/r, where r is the radius curvature
  // this is copied from ParticlePropagator
  // compute radius in [m]
  double radius = pt_in_gev / (charge * _bz) * 1.0E9/c_light;
  printf("radius: %f [m]\n", radius);
  // convert back to rho, in cm
  return -1 / (1e2 * radius);
}

double RaveConverter::getRhoAlt(double qoverp, double theta) {
  // double eta = -std::log(std::tan(theta/2));
  double cosh_eta = 1/std::sin(theta);
  double rho = - qoverp * _bz * cosh_eta * c_light * 1e-9 * 1e-2;
  return rho;
}
double RaveConverter::getQOverP(double rho, double theta) {
  return - rho * std::sin(theta) / (c_light * 1e-11 * _bz);
}
double RaveConverter::getDrhoDqoverp(double theta) {
  double cosh_eta = 1/std::sin(theta);
  return - _bz * cosh_eta * c_light * 1e-11;
}
double RaveConverter::getDrhoDtheta(double qoverp, double theta) {
  return c_light*1e-11 * _bz*qoverp / (std::tan(theta) * std::sin(theta));
}


std::vector<rave::Track> RaveConverter::getRaveTracks(
  const std::vector<Candidate*>& in) {
  std::vector<rave::Track> tracks;
  for (const auto& deltrack: in) {
    rave::Vector6D state = getState(deltrack);
    rave::PerigeeCovariance5D cov5d = getPerigeeCov(deltrack);
    int charge = deltrack->Charge;
    rave::Covariance6D cov6d = _to_rave.convert(cov5d, state, charge);
    rave::Track track(state, cov6d, charge, 0.0, 0.0, deltrack);
    tracks.push_back(track);
  }
  return tracks;
}


// define the utility functions which are forward declared above
namespace {
  // walk up the candidate tree to get the generated particle
  Candidate* get_part(Candidate* cand) {
    if (cand->GetCandidates()->GetEntriesFast() == 0) {
      return cand;
    }
    Candidate* mother = static_cast<Candidate*>(cand->GetCandidates()->At(0));
    return get_part(mother);
  }
  // hackattack (because we trust the shit outa non-const version)
  const Candidate* get_part(const Candidate* cand) {
    return get_part(const_cast<Candidate*>(cand));
  }
  void print_rave_track_info(const rave::Track& track) {
    std::cout << track << std::endl;
    auto perigee = track.perigeeParameters();
    std::cout << "ep: " << perigee.epsilon() << ", z: " << perigee.zp()
	      << ", rho: " << perigee.rho()
	      << ", theta: " << perigee.theta()
	      << ", phi: " << perigee.phip()
	      << std::endl;
  }
  void print_more_info(const Candidate* cand) {
    float d0 = cand->Dxy;
    assert(d0 == cand->trkPar[TrackParam::D0]);
    float z = cand->trkPar[TrackParam::Z0];
    std::cout << "d0: " << d0*0.1 << ", z: " << z*0.1
	      << ", qop: " << cand->trkPar[TrackParam::QOVERP]*10
	      << ", theta: " << cand->trkPar[TrackParam::THETA]
	      << ", phi: " << cand->trkPar[TrackParam::PHI] << std::endl;
    auto& mom = cand->Momentum;
    std::cout << "(" << cand->Xd*0.1 << ", " << cand->Yd*0.1 << ", "
	      << z*0.1 << " , " << mom.X() << ", " << mom.Y() << ", "
	      << mom.Z() << " )" << std::endl;
  }

  void print_track_info(const Candidate* cand) {
    const TLorentzVector& mom = cand->Momentum;
    const TLorentzVector& pos = cand->Position;
    const Candidate* mother = get_part(cand);
    const TLorentzVector& mpos = mother->Position;

    float d0 = cand->Dxy;
    assert(d0 == cand->trkPar[TrackParam::D0]);
    float phi = cand->trkPar[TrackParam::PHI];
    float z = cand->trkPar[TrackParam::Z0];
    float phi0 = phi - std::copysign(3.14159/2, 1);
    float x = d0 * std::cos(phi0);
    float y = d0 * std::sin(phi0);
    std::cout << "d0: " << d0 << ", z: " << z
	      << ", qop: " << cand->trkPar[TrackParam::QOVERP]
	      << ", theta: " << cand->trkPar[TrackParam::THETA]
	      << ", phi: " << cand->trkPar[TrackParam::PHI] << std::endl;
    std::cout << "dphi: " << (phi0 - std::atan2(cand->Yd, cand->Xd))/3.1415
	      << "pi " << std::endl;
    std::cout << "initial x, y, z, r: " << mpos.X() << " " << mpos.Y() << " "
	      << mpos.Z() << " " << mpos.Perp() << std::endl;
    std::cout << "det edge x, y, z, r: " << pos.X() << " " << pos.Y() << " "
	      << pos.Z() << " " << pos.Perp() << std::endl;
    std::cout << "used  x, y, z " << x << " " << y << " "
	      << cand->Zd << std::endl;
    std::cout << "d0: " << d0 << " mm, charge: " << cand->Charge << ", PID: "
	      << cand->PID << ", pt " << mom.Pt() << " GeV"<< std::endl;

    std::cout << "other x, y, z " << cand->Xd << " " << cand->Yd << " "
	      << cand->Zd << std::endl;
    std::cout << "momentum x, y, z " << mom.Px() << " " << mom.Py() << " "
	      << mom.Pz() << std::endl;
  }
  std::ostream& operator<<(std::ostream& os,
			   const rave::PerigeeParameters5D& pars) {
#define WRITE(PAR) os << #PAR << ": " << pars.PAR() << ", "
    WRITE(epsilon);
    WRITE(zp);
    WRITE(phip);
    WRITE(theta);
    WRITE(rho);
#undef WRITE
    return os;
  }
}

// end of RaveConverter
// ________________________________________________________________________

#else // if NO_RAVE is set

#include <iostream>

// dummy stand-in class
SecondaryVertexTagging::SecondaryVertexTagging() {}
SecondaryVertexTagging::~SecondaryVertexTagging(){}

void SecondaryVertexTagging::Init() {
  std::cerr << "** WARNING: Rave was not found! This is a dummy class that "
	    << "will run no secondary vertex tagging!" << std::endl;
  fJetInputArray = ImportArray(
    GetString("JetInputArray", "FastJetFinder/jets"));
  fItJetInputArray = fJetInputArray->MakeIterator();
}
void SecondaryVertexTagging::Process() {
  fItJetInputArray->Reset();
  Candidate* jet;
  while((jet = static_cast<Candidate*>(fItJetInputArray->Next())))
  {
    SecondaryVertex test;
    test.Lxy = -1;
    test.config = "zork";
    jet->secondaryVertices.push_back(test);
  }
}
void SecondaryVertexTagging::Finish() {}

#endif // check for NO_RAVE flag
