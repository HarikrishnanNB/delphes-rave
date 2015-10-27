#include "RaveConverter.hh"

#include "classes/DelphesClasses.h"


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
