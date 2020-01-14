//
// Created by pierre-yves on 14/01/20.
//

#ifndef MATHUTILS_INTEGRATION2D_H
#define MATHUTILS_INTEGRATION2D_H

#include <vector>
#include "Vector3d.h"

namespace mathutils {

  /**
  * Class for computing 2d integrations.
  */
  template <class Scalar>
  class Integration2d {

   public:

    /// Contructor of the class.
    Integration2d(Scalar (*F)(Scalar x), const int& nb_Gauss_points)
        : m_nb_Gauss_points(nb_Gauss_points){
      m_integrand = F;
    }

   protected:

    /// Function to be integrated
    Scalar (*m_integrand)(Scalar x) = nullptr;

    /// Number of Gauss points.
    int m_nb_Gauss_points;
    
  };

}

#endif //MATHUTILS_INTEGRATION2D_H
