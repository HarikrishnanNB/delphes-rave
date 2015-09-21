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


/** \class HDF5Writer
 *
 * writes jet information to HDF5
 *
 *  \author Dan Guest
 *
 */

#include "modules/HDF5Writer.h"
#include "external/h5/h5types.hh"
#include "external/flavortag/SecondaryVertex.hh"

#include "classes/DelphesClasses.h"
#include "ExRootAnalysis/ExRootConfReader.h"
// as a hack we get the output file name from ExRootTreeWriter
#include "ExRootAnalysis/ExRootTreeWriter.h"

#include "TObjArray.h"
#include "TFolder.h"

#include <iostream>
#include <string>
#include <limits>

namespace {
  // double inf = std::numeric_limits<double>::infinity();
  std::string remove_extension(const std::string& filename) {
    size_t lastdot = filename.find_last_of(".");
    if (lastdot == std::string::npos) return filename;
    return filename.substr(0, lastdot);
  }
}

//------------------------------------------------------------------------------

HDF5Writer::HDF5Writer() :
  fItInputArray(0), m_out_file(0), m_hl_jet_buffer(0),
  m_ml_jet_buffer(0)
{
}

//------------------------------------------------------------------------------

HDF5Writer::~HDF5Writer()
{
  delete m_out_file;
  delete m_hl_jet_buffer;
  delete m_ml_jet_buffer;
  delete fItInputArray;
}

//------------------------------------------------------------------------------
namespace out {
  // insering a compound type requires that `type(Class)` is defined
  H5::CompType type(Vertex);
  H5::CompType type(Jet);
}

void HDF5Writer::Init()
{
  fInputArray = ImportArray(
    GetString("JetInputArray", "UniqueObjectFinder/jets"));
  fItInputArray = fInputArray->MakeIterator();

  // get the name of the root output file
  auto* treeWriter = static_cast<ExRootTreeWriter*>(
    GetFolder()->FindObject("TreeWriter"));
  const std::string output_file = remove_extension(
    treeWriter->GetOutputFileName());

  // create the hdf5 output file
  std::string output_ext = GetString("OutputExtension", ".ntuple.h5");
  std::string hdf_out = output_file + output_ext;
  m_out_file = new H5::H5File(hdf_out, H5F_ACC_TRUNC);

  auto hl_jtype = out::type(out::HighLevelJet());
  auto ml_jtype = out::type(out::MediumLevelJet());

  m_hl_jet_buffer = new OneDimBuffer<out::HighLevelJet>(
    *m_out_file, "high_level_jets", hl_jtype, 1000);
  m_ml_jet_buffer = new OneDimBuffer<out::MediumLevelJet>(
    *m_out_file, "medium_level_jets", ml_jtype, 1000);

  // create the output text file
  std::string text_file_ext = GetString("TextFileExtension", "");
  if (text_file_ext.size() > 0) {
    m_output_stream.open(output_file + text_file_ext);
  }
}

namespace out {
  int simple_flavor(int flav) {
    switch (flav) {
    case 4: return 4;
    case 5: return 5;
    case 15: return 15;
    default: return 0;
    }
  }

  JetParameters::JetParameters(Candidate& jet):
    pt(jet.Momentum.Pt()),
    eta(jet.Momentum.Eta()),
    flavor(simple_flavor(jet.Flavor))
  {
  }
  HighLevelTracking::HighLevelTracking(const ::HighLevelTracking& hlTrk):
    track_2_d0_significance(hlTrk.track2d0sig),
    track_3_d0_significance(hlTrk.track3d0sig),
    track_2_z0_significance(hlTrk.track2z0sig),
    track_3_z0_significance(hlTrk.track3z0sig),
    n_tracks_over_d0_threshold(hlTrk.tracksOverIpThreshold),
    jet_prob(hlTrk.jetProb),
    jet_width_eta(hlTrk.jetWidthEta),
    jet_width_phi(hlTrk.jetWidthPhi)
  {
  }
  HighLevelSecondaryVertex::HighLevelSecondaryVertex(
    const ::HighLevelSvx& hlSvx):
    vertex_significance(hlSvx.Lsig),
    n_secondary_vertices(hlSvx.NVertex),
    n_secondary_vertex_tracks(hlSvx.NTracks),
    delta_r_vertex(hlSvx.DrJet),
    vertex_mass(hlSvx.Mass),
    vertex_energy_fraction(hlSvx.EnergyFraction)
  {
  }

  HighLevelJet::HighLevelJet(Candidate& jet):
    jet_parameters(jet),
    tracking(jet.hlTrk),
    vertex(jet.hlSvx)
  {
  }

  VertexTrack::VertexTrack(const SecondaryVertexTrack& tk):
    d0(tk.d0), z0(tk.z0),
    d0_uncertainty(tk.d0err), z0_uncertainty(tk.z0err),
    pt(tk.pt),
    delta_phi_jet(tk.dphi), delta_eta_jet(tk.deta)
  {
  }
  SecondaryVertex::SecondaryVertex(const ::SecondaryVertex& vx):
    mass(vx.mass),
    displacement(vx.Mag()),
    delta_eta_jet(vx.deta),
    delta_phi_jet(vx.dphi),
    displacement_significance(vx.Lsig)
  {
    for (const auto& track: vx.tracks_along_jet) {
      associated_tracks.push_back(track);
    }
  }
  MediumLevelJet::MediumLevelJet(Candidate& jet):
    jet_parameters(jet)
  {
    for (const auto& trk: jet.primaryVertexTracks) {
      primary_vertex_tracks.push_back(trk);
    }
    for (const auto& vx: jet.secondaryVertices) {
      secondary_vertices.push_back(SecondaryVertex(vx));
    }
  }

  // ____________________________________________________________________
  // HDF5 types

  // insering a compound type requires that `type(Class)` is defined
  // high level variables
  H5::CompType type(JetParameters) {
    H5::CompType out(sizeof(JetParameters));
    H5_INSERT(out, JetParameters, pt);
    H5_INSERT(out, JetParameters, eta);
    H5_INSERT(out, JetParameters, flavor);
    return out;
  }

  H5::CompType type(HighLevelTracking) {
    H5::CompType out(sizeof(HighLevelTracking));
    H5_INSERT(out, HighLevelTracking, track_2_d0_significance);
    H5_INSERT(out, HighLevelTracking, track_3_d0_significance);
    H5_INSERT(out, HighLevelTracking, track_2_z0_significance);
    H5_INSERT(out, HighLevelTracking, track_3_z0_significance);
    H5_INSERT(out, HighLevelTracking, n_tracks_over_d0_threshold);
    H5_INSERT(out, HighLevelTracking, jet_prob);
    H5_INSERT(out, HighLevelTracking, jet_width_eta);
    H5_INSERT(out, HighLevelTracking, jet_width_phi);
    return out;
  }

  H5::CompType type(HighLevelSecondaryVertex) {
    H5::CompType out(sizeof(HighLevelSecondaryVertex));
    H5_INSERT(out, HighLevelSecondaryVertex, vertex_significance);
    H5_INSERT(out, HighLevelSecondaryVertex, n_secondary_vertices);
    H5_INSERT(out, HighLevelSecondaryVertex, n_secondary_vertex_tracks);
    H5_INSERT(out, HighLevelSecondaryVertex, delta_r_vertex);
    H5_INSERT(out, HighLevelSecondaryVertex, vertex_mass);
    H5_INSERT(out, HighLevelSecondaryVertex, vertex_energy_fraction);
    return out;
  }
  H5::CompType type(HighLevelJet) {
    H5::CompType out(sizeof(HighLevelJet));
    H5_INSERT(out, HighLevelJet, jet_parameters);
    H5_INSERT(out, HighLevelJet, tracking);
    H5_INSERT(out, HighLevelJet, vertex);
    return out;
  }

  // medium-level variables
  H5::CompType type(VertexTrack) {
    H5::CompType out(sizeof(VertexTrack));
    H5_INSERT(out, VertexTrack, d0);
    H5_INSERT(out, VertexTrack, z0);
    H5_INSERT(out, VertexTrack, d0_uncertainty);
    H5_INSERT(out, VertexTrack, z0_uncertainty);
    H5_INSERT(out, VertexTrack, pt);
    H5_INSERT(out, VertexTrack, delta_eta_jet);
    H5_INSERT(out, VertexTrack, delta_phi_jet);
    return out;
  }
  H5::CompType type(SecondaryVertex) {
    H5::CompType out(sizeof(SecondaryVertex));
    H5_INSERT(out, SecondaryVertex, mass);
    H5_INSERT(out, SecondaryVertex, displacement);
    H5_INSERT(out, SecondaryVertex, delta_eta_jet);
    H5_INSERT(out, SecondaryVertex, delta_phi_jet);
    H5_INSERT(out, SecondaryVertex, displacement_significance);
    H5_INSERT(out, SecondaryVertex, associated_tracks);
    return out;
  }
  H5::CompType type(MediumLevelJet) {
    H5::CompType out(sizeof(MediumLevelJet));
    H5_INSERT(out, MediumLevelJet, jet_parameters);
    H5_INSERT(out, MediumLevelJet, primary_vertex_tracks);
    H5_INSERT(out, MediumLevelJet, secondary_vertices);
    return out;
  }
}

//------------------------------------------------------------------------------

void HDF5Writer::Finish()
{
  m_hl_jet_buffer->flush();
  m_hl_jet_buffer->close();

  m_ml_jet_buffer->flush();
  m_ml_jet_buffer->close();

  if (m_output_stream.is_open()) {
    m_output_stream.close();
  }
}

//------------------------------------------------------------------------------

void HDF5Writer::Process()
{
  fItInputArray->Reset();
  Candidate* jet;
  while ((jet = static_cast<Candidate*>(fItInputArray->Next()))) {
    out::MediumLevelJet ml_jet(*jet);
    if (m_output_stream.is_open()) {
      m_output_stream << ml_jet << "\n";
    }
    m_hl_jet_buffer->push_back(*jet);
    m_ml_jet_buffer->push_back(ml_jet);
  }
}

//------------------------------------------------------------------------------
//
namespace out {
// ostream operators
  std::ostream& operator<<(std::ostream& out, const JetParameters& pars) {
    out << pars.pt << ", ";
    out << pars.eta << ", ";
    out << pars.flavor;
    return out;
  }
  std::ostream& operator<<(std::ostream& out, const HighLevelTracking& pars) {
    out << pars.track_2_d0_significance << ", ";
    out << pars.track_3_d0_significance << ", ";
    out << pars.track_2_z0_significance << ", ";
    out << pars.track_3_z0_significance << ", ";
    out << pars.n_tracks_over_d0_threshold << ", ";
    out << pars.jet_prob << ", ";
    out << pars.jet_width_eta << ", ";
    out << pars.jet_width_phi;
    return out;
  }
  std::ostream& operator<<(std::ostream& out,
			   const HighLevelSecondaryVertex& pars) {
    out << pars.vertex_significance << ", ";
    out << pars.n_secondary_vertices << ", ";
    out << pars.n_secondary_vertex_tracks << ", ";
    out << pars.delta_r_vertex << ", ";
    out << pars.vertex_mass << ", ";
    out << pars.vertex_energy_fraction;
    return out;
  }
  std::ostream& operator<<(std::ostream& out, const HighLevelJet& pars) {
    out << pars.jet_parameters << ", " << pars.tracking << ", "
	<< pars.vertex;
    return out;
  }
  std::ostream& operator<<(std::ostream& out, const VertexTrack& pars) {
    out << pars.d0 << ", ";
    out << pars.z0 << ", ";
    out << pars.d0_uncertainty << ", ";
    out << pars.z0_uncertainty << ", ";
    out << pars.pt << ", ";
    out << pars.delta_phi_jet << ", ";
    out << pars.delta_eta_jet;
    return out;
  }
  std::ostream& operator<<(std::ostream& out, const SecondaryVertex& pars) {
    out << pars.mass << ", ";
    out << pars.displacement << ", ";
    out << pars.delta_eta_jet << ", ";
    out << pars.delta_phi_jet << ", ";
    out << pars.displacement_significance << ", [";
    size_t n_trk = pars.associated_tracks.size();
    for (size_t iii = 0; iii < n_trk; iii++) {
      const auto& trk = pars.associated_tracks.at(iii);
      out << "{" << trk << "}";
      if (iii != (n_trk - 1)) out << ", ";
    }
    out << "]";
    return out;
  }
  std::ostream& operator<<(std::ostream& out, const MediumLevelJet& pars) {
    out << pars.jet_parameters << ", [";
    size_t n_trk = pars.primary_vertex_tracks.size();
    for (size_t iii = 0; iii < n_trk; iii++) {
      const auto& trk = pars.primary_vertex_tracks.at(iii);
      out << "{" << trk << "}";
      if (iii != (n_trk - 1)) out << ", ";
    }
    out << "], [";
    size_t n_svx = pars.secondary_vertices.size();
    for (size_t iii = 0; iii < n_svx; iii++) {
      const auto& vx = pars.secondary_vertices.at(iii);
      out << "{" << vx << "}";
      if (iii != (n_svx - 1)) out << ", ";
    }
    out << "]";
    return out;
  }
}
