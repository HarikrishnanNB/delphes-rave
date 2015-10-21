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

#ifndef IPCovSmearing_h
#define IPCovSmearing_h

/** \class IPCovSmearing
 *
 *  Performs transverse impact parameter smearing.
 *
 *  \author M. Selvaggi - UCL, Louvain-la-Neuve
 *
 */

#include "classes/DelphesModule.h"
#include <vector>

#ifndef __CINT__
#include <unordered_map>
#include <random>
#include <Eigen/Cholesky>
typedef Eigen::Matrix<double, 5, 5> CovMatrix;
typedef Eigen::Matrix<double, 5, 1> TrackVector;
#endif

class TIterator;
class TObjArray;
class DelphesFormula;
class TFile;


class IPCovSmearing: public DelphesModule
{
public:

  IPCovSmearing();
  ~IPCovSmearing();

  void Init();
  void Process();
  void Finish();

private:

  std::vector<double> ptbins, etabins;
  unsigned long long fNBinMisses;

#ifndef __CINT__
  std::default_random_engine fRandomGenerator;

  // unility for bins that aren't covered for our smearing
  // this is to grab one that is
  std::pair<int,int> getValidBins(int ptbin, int etabin);
  TrackVector getRandomVector();
  std::unordered_map<
    int, std::unordered_map<int, CovMatrix> > fSmearingMatrices;
  std::unordered_map<
    int, std::unordered_map<int, CovMatrix> > fCovarianceMatrices;
#endif

  TIterator *fItInputArray; //!

  const TObjArray *fInputArray; //!

  TObjArray *fOutputArray; //!

  ClassDef(IPCovSmearing, 1)
};

#endif
