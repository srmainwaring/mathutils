//
// Created by pierre-yves on 14/01/20.
//

#ifndef MATHUTILS_INTEGRATION2D_H
#define MATHUTILS_INTEGRATION2D_H

#include "QuadratureTables.h"

namespace mathutils {

  /**
  * Class for computing 2d integrations based on a Gauss quadrature.
  */
  template<typename Tin, typename Tout>
  class Integration2d {

   public:

    /// Contructor of the class.
    Integration2d(Tout (*F)(Tin x), const int& order)
        : m_order(order){
      m_integrand = F;
      m_tables = QuadratureTables();
    }

    /// This function computes the 2d integration.
    virtual Tout Compute() = 0;

   protected:

    /// Function to be integrated
    Tout (*m_integrand)(Tin x) = nullptr;

    /// Order of the quadrature.
    int m_order;

    /// Table of quadrature coefficients.
    QuadratureTables m_tables;

  };

}

#endif //MATHUTILS_INTEGRATION2D_H
