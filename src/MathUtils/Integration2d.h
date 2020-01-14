//
// Created by pierre-yves on 14/01/20.
//

#ifndef MATHUTILS_INTEGRATION2D_H
#define MATHUTILS_INTEGRATION2D_H

#include <vector>

namespace mathutils {

  /**
  * Class for computing 2d integrations.
  */
  template<typename Tin, typename Tout>
  class Integration2d {

   public:

    /// Contructor of the class.
    Integration2d(Tout (*F)(Tin x), const int& nb_Gauss_points)
        : m_nb_Gauss_points(nb_Gauss_points){
      m_integrand = F;
    }

   protected:

    /// Function to be integrated
    Tout (*m_integrand)(Tin x) = nullptr;

    /// Number of Gauss points.
    int m_nb_Gauss_points;

  };

}

#endif //MATHUTILS_INTEGRATION2D_H
