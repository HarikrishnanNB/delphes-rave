#include "hl_vars.hh"

// #include "classes/DelphesClasses.h"
#include "constants_jetprob.hh"
#include "enums_track.hh"
#include "math.hh"

#include <vector>
#include <utility>
#include <cassert>
#include <cmath>
#include <algorithm>
#include <limits>

// #include <iostream>

struct TrackParameters;

namespace {
  typedef std::vector<TrackParameters> Tracks;

  const double pi = std::atan2(0, -1);
  static_assert(std::numeric_limits<double>::has_infinity, "need inf");
  const double inf = std::numeric_limits<double>::infinity();
  static_assert(std::numeric_limits<double>::has_quiet_NaN, "need NaN");
  const double NaN = std::numeric_limits<double>::quiet_NaN();
  template<typename T>
  bool by_descending_first(std::pair<double, T> v1, std::pair<double, T> v2) {
    return v1.first > v2.first;
  }
  double get_jet_prob(const std::vector<TrackParameters>&);

  // see hardcoded parameters in constants_jetprob.hh
  double get_track_prob(double d0sig);

  typedef std::pair<double, double> JetWidth;
  JetWidth jet_width2_eta_phi(const TVector3& jet, const Tracks& tracks);

}

// __________________________________________________________________________
// SVX


HighLevelSvx::HighLevelSvx():
  Lsig(NaN), NVertex(-1), NTracks(-1), DrJet(NaN), Mass(NaN),
  EnergyFraction(NaN)
{
}

void HighLevelSvx::fill(const TVector3& jvec,
			const std::vector<SecondaryVertex>& vertices,
			size_t skip) {
  double sum_dr_tracks = 0;
  int sum_tracks = 0;
  int sum_vertices = 0;
  double sum_mass = 0;
  double sum_efrac = 0;

  // copied these variables from jetfitter
  double sum_sig = 0;
  double sum_inverr = 0;
  const size_t n_vtx = vertices.size();
  for (size_t vxn = skip; vxn < n_vtx; vxn++) {
    const auto& vx = vertices.at(vxn);
    if (vx.Pt() == 0) continue;
    sum_vertices++;
    double delta_r = jvec.DeltaR(vx);
    sum_dr_tracks += vx.nTracks * delta_r;
    sum_tracks += vx.nTracks;

    sum_mass += vx.mass;
    sum_efrac += vx.eFrac;

    sum_sig += vx.Mag() / vx.decayLengthVariance;
    sum_inverr += 1/vx.decayLengthVariance;
  }
  bool has_vx = (sum_vertices > 0) && (sum_tracks > 0);

  // save summary info
  Lsig = has_vx ? sum_sig / sqrt(sum_inverr) : -1;
  NVertex = has_vx ? sum_vertices               : -1;
  NTracks = has_vx ? sum_tracks                 : -1;
  DrJet = has_vx ? sum_dr_tracks / sum_tracks : inf;
  Mass = has_vx ? sum_mass                   : -1;
  EnergyFraction = has_vx ? sum_efrac : -inf;
}

std::ostream& operator<<(std::ostream& os, const HighLevelSvx& hl) {
#define DUMP(VAR) os << #VAR ": " << hl.VAR << ", "
  DUMP(Lsig);
  DUMP(NVertex);
  DUMP(NTracks);
  DUMP(DrJet);
  DUMP(Mass);
  DUMP(EnergyFraction);
#undef DUMP
  return os;
}

// ________________________________________________________________________
// Tracking

TrackParameters::TrackParameters(const float trkPar[5],
				 const float trkCov[15]):
  d0(trkPar[trk::D0]),
  z0(trkPar[trk::Z0]),
  phi(trkPar[trk::PHI]),
  theta(trkPar[trk::THETA]),
  qoverp(trkPar[trk::QOVERP]),
  d0err(std::sqrt(trkCov[trk::D0D0])),
  z0err(std::sqrt(trkCov[trk::Z0Z0]))
{
}

std::ostream& operator<<(std::ostream& os, const TrackParameters& hl) {
#define DUMP(VAR) os << #VAR ": " << hl.VAR << ", "
  DUMP(d0);
  DUMP(z0);
  DUMP(phi);
  DUMP(theta);
  DUMP(d0err);
  DUMP(z0err);
#undef DUMP
  return os;
}

HighLevelTracking::HighLevelTracking():
  track2d0sig(NaN), track3d0sig(NaN),
  track2z0sig(NaN), track3z0sig(NaN),
  tracksOverIpThreshold(-1), jetProb(NaN),
  jetWidthEta(NaN), jetWidthPhi(NaN)
{
}

void HighLevelTracking::fill(const TVector3& jet,
			     const std::vector<TrackParameters>& pars,
			     double ip_threshold) {
  // size_t n_tracks = pars.size();
  double jet_phi = jet.Phi();
  assert(std::abs(jet_phi) <= pi);
  std::vector<std::pair<double, TrackParameters> > tracks_by_ip;

  // zero some things
  track2d0sig = -inf;
  track2z0sig = -inf;
  track3d0sig = -inf;
  track3z0sig = -inf;
  tracksOverIpThreshold = 0;
  jetProb = -1;
  jetWidthEta = -inf;
  jetWidthPhi = -inf;

  auto eta_phi = jet_width2_eta_phi(jet, pars);
  if (eta_phi.first > 0) {
    jetWidthEta = std::sqrt(eta_phi.first);
    jetWidthPhi = std::sqrt(eta_phi.second);
  }

  if (pars.size() == 0) return;
  jetProb = get_jet_prob(pars);

  // what follows uses numbered tracks (track counting)
  if (pars.size() < 2) return;

  for (const auto& par: pars) {
    double diff = std::abs(jet_phi - par.phi);
    int sign = (diff > 3*pi/4 || diff < pi/2) ? 1 : -1;
    double ip = std::copysign(par.d0, sign);
    tracks_by_ip.emplace_back(ip, par);
    if ((ip / par.d0err) > ip_threshold) tracksOverIpThreshold++;
  }

  std::sort(tracks_by_ip.begin(), tracks_by_ip.end(),
	    by_descending_first<TrackParameters>);
  {
    const auto& trk2 = tracks_by_ip.at(1);
    const auto& p2 = trk2.second;
    track2d0sig = std::copysign(p2.d0 / p2.d0err, trk2.first);
    track2z0sig = std::abs(p2.z0 / p2.z0err);
  }
  if (tracks_by_ip.size() < 3) return;
  {
    const auto& trk3 = tracks_by_ip.at(2);
    const auto& p3 = trk3.second;
    track3d0sig = std::copysign(p3.d0 / p3.d0err, trk3.first);
    track3z0sig = std::abs(p3.z0 / p3.z0err);
  }
}

std::ostream& operator<<(std::ostream& os, const HighLevelTracking& hl) {
#define DUMP(VAR) os << #VAR ": " << hl.VAR << ", "
  DUMP(track2d0sig);
  DUMP(track2z0sig);
  DUMP(track3d0sig);
  DUMP(track3z0sig);
  DUMP(tracksOverIpThreshold);
  DUMP(jetProb);
  DUMP(jetWidthEta);
  DUMP(jetWidthPhi);
#undef DUMP
  return os;
}

namespace {
  double gauss_prob(double sig, const double norm, const double width) {
    const double sqp = std::sqrt(pi);
    const double sq2 = std::sqrt(2);
    return sqp / 2 * norm * width * (1 - std::erf( sig / (sq2 * width) ));
  }
  double exp_prob(double sig, const double off, const double mult) {
    return 1 / mult * std::exp(-off - mult*sig);
  }
  double get_track_prob(double sig) {
    using namespace jetprob;
    double prob = gauss_prob(sig, P0, P1) + gauss_prob(sig, P2, P3) +
      exp_prob(sig, P4, P5) + exp_prob(sig, P6, P7);
    return prob;
  }
  double get_jet_prob(const std::vector<TrackParameters>& pars) {
    double p0 = 1.0;
    for (const auto& par: pars) {
      double sig = par.d0 / par.d0err;
      double prob = get_track_prob(std::abs(sig));
      p0 *= prob;
    }
    int n_trk = pars.size();
    double corrections = 0;
    for (int k = 0; k < n_trk; k++) {
      corrections += std::pow( -std::log(p0), k) / std::tgamma(k + 1);
    }
    return p0 * corrections;
  }
  JetWidth jet_width2_eta_phi(const TVector3& jet, const Tracks& tracks) {
    if (tracks.size() < 1) return {-1, -1};

    const double jet_eta = jet.Eta();
    const double jet_phi = jet.Phi();
    double sum_pt_times_eta2 = 0;
    double sum_pt_times_phi2 = 0;
    double sum_pt = 0;
    for (const auto& trk: tracks) {
      double eta = -std::log(std::tan(trk.theta/2));
      double deta = eta - jet_eta;
      double dphi = phi_mpi_pi(trk.phi, jet_phi);
      double track_pt = std::abs(1 / (trk.qoverp * std::cosh(eta)));
      sum_pt += track_pt;
      sum_pt_times_eta2 += track_pt * deta*deta;
      sum_pt_times_phi2 += track_pt * dphi*dphi;
    }
    assert(sum_pt > 0);
    return {sum_pt_times_eta2 / sum_pt, sum_pt_times_phi2 / sum_pt};
  }
}
