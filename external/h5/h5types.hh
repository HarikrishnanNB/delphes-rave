#ifndef H5_TYPES_HH
#define H5_TYPES_HH

#include "h5container.hh"

#include "H5Cpp.h"
#include <vector>
#include <string>

#define H5_INSERT(INTO, CLASS, MEMBER)					\
  h5::insert(INTO, #MEMBER, offsetof(CLASS, MEMBER), &CLASS::MEMBER)


namespace h5 {
  H5::DataType type(int);
  H5::DataType type(double);
  H5::DataType type(float);
  template <typename T>
  H5::DataType type(const h5::vector<T>&) {
    const auto subtype = type(T());
    return H5::VarLenType(&subtype);
  }

  template <typename M, typename T>
  void insert(H5::CompType& into, const std::string& name,
	      size_t offset, M T::*) {
    into.insertMember(name, offset, type(M()));
  }

}

#endif
