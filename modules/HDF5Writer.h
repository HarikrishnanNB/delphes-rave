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

#ifndef HDF5Writer_h
#define HDF5Writer_h

/** \class HDF5Writer
 *
 *  Selects candidates from the InputArray according to the efficiency formula.
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */



class TObjArray;
class DelphesFormula;
class Candidate;
class SecondaryVertexTrack;
class SecondaryVertex;
class HighLevelTracking;
class HighLevelSvx;

#include <fstream>

#ifndef __CINT__

#include "external/h5/OneDimBuffer.hh"
#include "external/h5/h5container.hh"

#include "classes/DelphesModule.h"
#include "h5/bork.hh"

#include "H5Cpp.h"

namespace out {

  // high-level
  struct JetParameters {
    JetParameters(Candidate& jet);
    JetParameters() = default;
    double pt;
    double eta;
    int flavor;
  };
  H5::CompType type(JetParameters);
  std::ostream& operator<<(std::ostream&, const JetParameters&);

  struct HighLevelTracking {
    HighLevelTracking(const ::HighLevelTracking&);
    HighLevelTracking() = default;
    double track_2_d0_significance;
    double track_3_d0_significance;
    double track_2_z0_significance;
    double track_3_z0_significance;
    int    n_tracks_over_d0_threshold;
    double jet_prob;
    double jet_width_eta;
    double jet_width_phi;
  };
  H5::CompType type(HighLevelTracking);
  std::ostream& operator<<(std::ostream&, const HighLevelTracking&);

  struct HighLevelSecondaryVertex {
    HighLevelSecondaryVertex(const ::HighLevelSvx&);
    HighLevelSecondaryVertex() = default;
    double vertex_significance;
    int    n_secondary_vertices;
    int    n_secondary_vertex_tracks;
    double delta_r_vertex;
    double vertex_mass;
    double vertex_energy_fraction;
  };
  H5::CompType type(HighLevelSecondaryVertex);
  std::ostream& operator<<(std::ostream&, const HighLevelSecondaryVertex&);

  struct HighLevelJet {
    HighLevelJet(Candidate& jet);
    HighLevelJet() = default;
    // basic parameters
    JetParameters jet_parameters;

    // track-based
    HighLevelTracking tracking;

    // secondary vertex
    HighLevelSecondaryVertex vertex;
  };
  H5::CompType type(HighLevelJet);
  std::ostream& operator<<(std::ostream&, const HighLevelJet&);

  // medium level
  struct VertexTrack {
    VertexTrack(const SecondaryVertexTrack&);
    VertexTrack() = default;
    double d0;
    double z0;
    double d0_uncertainty;
    double z0_uncertainty;
    double pt;
    double delta_phi_jet;
    double delta_eta_jet;
    double weight;
  };
  H5::CompType type(VertexTrack);
  std::ostream& operator<<(std::ostream&, const VertexTrack&);
  std::ostream& operator<<(std::ostream&, const h5::vector<VertexTrack>&);

  struct SecondaryVertex {
    SecondaryVertex(const ::SecondaryVertex&);
    SecondaryVertex() = default;
    double mass;
    double displacement;
    double delta_eta_jet;
    double delta_phi_jet;
    double displacement_significance;
    h5::vector<VertexTrack> associated_tracks;
  };
  H5::CompType type(SecondaryVertex);
  std::ostream& operator<<(std::ostream&, const SecondaryVertex&);
  std::ostream& operator<<(std::ostream&, const h5::vector<SecondaryVertex>&);

  struct MediumLevelJet {
    MediumLevelJet(Candidate& jet);
    MediumLevelJet() = default;
    JetParameters jet_parameters;

    h5::vector<VertexTrack> primary_vertex_tracks;
    h5::vector<SecondaryVertex> secondary_vertices;
  };
  H5::CompType type(MediumLevelJet);
  std::ostream& operator<<(std::ostream&, const MediumLevelJet&);

  // this one is only used for the ostream
  struct SuperJet {
    SuperJet(Candidate& jet);
    SuperJet() = default;
    JetParameters jet_parameters;

    // high-level
    HighLevelTracking tracking;
    HighLevelSecondaryVertex vertex;
    // medium level
    h5::vector<VertexTrack> primary_vertex_tracks;
    h5::vector<SecondaryVertex> secondary_vertices;
  };
  std::ostream& operator<<(std::ostream&, const SuperJet&);
}

#else  // CINT include dummy

namespace H5 {
  class H5File;
}

#endif

class HDF5Writer: public DelphesModule
{
public:

  HDF5Writer();
  ~HDF5Writer();

  void Init();
  void Process();
  void Finish();

private:

  TIterator *fItInputArray; //!

  const TObjArray *fInputArray; //!

  H5::H5File* m_out_file;

  double fPTMin;
  double fAbsEtaMax;

#ifndef __CINT__
  OneDimBuffer<out::HighLevelJet>* m_hl_jet_buffer;
  OneDimBuffer<out::MediumLevelJet>* m_ml_jet_buffer;
#endif
  std::ofstream m_output_stream;

  ClassDef(HDF5Writer, 1)
};

#endif
