#ifndef RAVE_CONVERTER_HH
#define RAVE_CONVERTER_HH

#include "rave/Vector6D.h"
#include "rave/Track.h"

#include "RaveBase/Converters/interface/PerigeeToRaveObjects.h"
#include "RaveBase/Converters/interface/RaveToPerigeeObjects.h"

#include <vector>

class Candidate;

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

#endif
