#include "h5types.hh"

namespace h5 {
  H5::DataType type(int) { return H5::PredType::NATIVE_INT; }
  H5::DataType type(double) {return H5::PredType::NATIVE_DOUBLE; }
  H5::DataType type(float) {return H5::PredType::NATIVE_FLOAT; }
}
