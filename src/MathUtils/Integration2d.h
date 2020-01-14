//
// Created by pierre-yves on 14/01/20.
//

#ifndef MATHUTILS_INTEGRATION2D_H
#define MATHUTILS_INTEGRATION2D_H

#include "QuadratureTables.h"
#include <memory>
#include "Vector3d.h"

namespace mathutils {

  /**
  * Class for computing 2d integrations based on a Gauss quadrature.
  */
  template<typename Tout>
  class Integration2d {

   public:

    /// Contructor of the class.
    Integration2d(Tout (*F)(Vector3d<double> x), const int& order)
        : m_order(order){
      m_integrand = F;
      m_tables = std::make_unique<QuadratureTables>();
    }

    /// This function computes the 2d integration.
    virtual Tout Compute(const Vector3d<double> &vertex_1, const Vector3d<double> &vertex_2,
                         const Vector3d<double> &vertex_3) = 0;

   protected:

    /// Function to be integrated
    Tout (*m_integrand)(Vector3d<double> x) = nullptr;

    /// Order of the quadrature.
    int m_order;

    /// Table of quadrature coefficients.
    std::unique_ptr<QuadratureTables> m_tables;

  };

}

#endif //MATHUTILS_INTEGRATION2D_H
