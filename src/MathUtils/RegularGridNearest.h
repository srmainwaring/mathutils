//
// Created by frongere on 18/05/22.
//

#ifndef MATHUTILS_REGULARGRIDNEAREST_H
#define MATHUTILS_REGULARGRIDNEAREST_H

#include <vector>
#include <array>
#include <memory>
#include <iostream>
#include <boost/multi_array.hpp>

namespace mathutils {

  template<typename Data_t, size_t _dim, typename Coord_t=double>
  struct GridNode {
   public:
    GridNode(const std::array<Coord_t, _dim> &coords, const Data_t &val) :
        m_coords(coords),
        m_val(val) {}

    std::array<Coord_t, _dim> coords() const { return m_coords; }

    Data_t val() const { return m_val; }

   private:
    std::array<Coord_t, _dim> m_coords;
    Data_t m_val;

  };


  // TODO: beaucoup de recoupement a factoriser avec RegularGridInterpolation -> classe de base a mettre en place

  template<typename Data_t, size_t _dim, typename Coord_t=double>
  class RegularGridNearest {

    using GridType = std::vector<std::vector<Coord_t>>;
    using NDArray = boost::multi_array<Data_t, _dim>;

   public:

    RegularGridNearest();

    RegularGridNearest(const RegularGridNearest &other);

    RegularGridNearest &operator=(const RegularGridNearest &other);

    void AddCoord(const std::vector<Coord_t> &vec);

    void AddVar(const NDArray &ndarray);

    std::vector<Coord_t> GetCoord(unsigned int idx) const;

    Coord_t GetCoordMinVal(unsigned int idx) const;

    Coord_t GetCoordMaxVal(unsigned int idx) const;

    size_t GetDimension() const;

    GridNode<Data_t, _dim, Coord_t> Nearest(const std::array<Coord_t, _dim> &point) const;

    bool IsValid() const;

   private:
    GridType m_grid_list;
    NDArray m_data;

  };

  template<typename Data_t, size_t _dim, typename Coord_t>
  RegularGridNearest<Data_t, _dim, Coord_t>::RegularGridNearest() {
    m_grid_list.reserve(_dim);
  }

  template<typename Data_t, size_t _dim, typename Coord_t>
  RegularGridNearest<Data_t, _dim, Coord_t>::RegularGridNearest(const RegularGridNearest &other) {
    *this = other;
  }

  template<typename Data_t, size_t _dim, typename Coord_t>
  RegularGridNearest<Data_t, _dim, Coord_t> &
  RegularGridNearest<Data_t, _dim, Coord_t>::operator=(const RegularGridNearest &other) {
    m_grid_list = other.m_grid_list;

    std::vector<size_t> ex;
    const size_t *shape = other.m_data.shape();
    ex.assign(shape, shape + other.m_data.num_dimensions());
    m_data.resize(ex);
    m_data = other.m_data;

    return *this;
  }

  template<typename Data_t, size_t _dim, typename Coord_t>
  void RegularGridNearest<Data_t, _dim, Coord_t>::AddCoord(const std::vector<Coord_t> &vec) {
    if (m_grid_list.size() == _dim) {
      std::cerr << "Attempting to add too much dimension to the interpolator" << std::endl;
      exit(EXIT_FAILURE);
    }
    m_grid_list.push_back(vec);
  }

  template<typename Data_t, size_t _dim, typename Coord_t>
  void RegularGridNearest<Data_t, _dim, Coord_t>::AddVar(const RegularGridNearest::NDArray &ndarray) {
// Verifying that every dimension has been filled
    if (m_grid_list.size() != _dim) {
      std::cerr << "Grid is incompletely defined" << std::endl;
      exit(EXIT_FAILURE);
    }

    // Copying the multi_array
    std::vector<size_t> ex;
    const size_t *shape = ndarray.shape();
    ex.assign(shape, shape + ndarray.num_dimensions());
    m_data.resize(ex);
    m_data = ndarray;
  }

  template<typename Data_t, size_t _dim, typename Coord_t>
  std::vector<Coord_t> RegularGridNearest<Data_t, _dim, Coord_t>::GetCoord(unsigned int idx) const {
    if (idx >= _dim) {
      std::cerr << "In RegularGridInterpolator, requesting a coordinate with index that is out of dimension"
                << std::endl;
      exit(EXIT_FAILURE);
    }
    return m_grid_list[idx];
  }

  template<typename Data_t, size_t _dim, typename Coord_t>
  Coord_t RegularGridNearest<Data_t, _dim, Coord_t>::GetCoordMinVal(unsigned int idx) const {
    if (idx >= _dim) {
      std::cerr << "In RegularGridInterpolator, requesting a coordinate with index that is out of dimension"
                << std::endl;
      exit(EXIT_FAILURE);
    }
    return *std::min_element(std::begin(m_grid_list[idx]), std::end(m_grid_list[idx]));
  }

  template<typename Data_t, size_t _dim, typename Coord_t>
  Coord_t RegularGridNearest<Data_t, _dim, Coord_t>::GetCoordMaxVal(unsigned int idx) const {
    if (idx >= _dim) {
      std::cerr << "In RegularGridInterpolator, requesting a coordinate with index that is out of dimension"
                << std::endl;
      exit(EXIT_FAILURE);
    }
    return *std::max_element(std::begin(m_grid_list[idx]), std::end(m_grid_list[idx]));
  }

  template<typename Data_t, size_t _dim, typename Coord_t>
  size_t RegularGridNearest<Data_t, _dim, Coord_t>::GetDimension() const {
    return _dim;
  }

  template<typename Data_t, size_t _dim, typename Coord_t>
  GridNode<Data_t, _dim, Coord_t>
  RegularGridNearest<Data_t, _dim, Coord_t>::Nearest(const std::array<Coord_t, _dim> &point) const {
    // TODO !!
//    Data_t val;

    std::array<unsigned int, _dim> indices;
    std::array<Coord_t, _dim> coords;

    for (unsigned int i = 0; i < _dim; ++i) {
      auto xi = point[i];

      int icell;

      // Looking for the cell in which the xi point lies
      const auto &coord = m_grid_list[i];
      if (xi < *coord.begin()) {
        icell = 0;
        // TODO
      } else if (xi >= (*(coord.end() - 1))) {
        icell = coord.size() - 1;
        // TODO
      } else {
        auto i_upper = std::upper_bound(coord.begin(), coord.end(), xi);
        icell = i_upper - coord.begin() - 1;

        if (xi >= 0.5 * (coord[icell + 1] + coord[icell])) {
          ++icell;
        }

        // TODO: tester si on prend icell ou icell+1 suivant la proximite de xi

      }
      indices[i] = icell;
      coords[i] = coord[icell];
    }

    Data_t val = m_data(indices);
    return GridNode<Data_t, _dim, Coord_t>(coords, val);
//    return m_data(indices);
  }

  template<typename Data_t, size_t _dim, typename Coord_t>
  bool RegularGridNearest<Data_t, _dim, Coord_t>::IsValid() const {
    return true;
  }

}  // mathutils


#endif //MATHUTILS_REGULARGRIDNEAREST_H
