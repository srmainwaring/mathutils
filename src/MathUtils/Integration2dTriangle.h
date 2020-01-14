//
// Created by pierre-yves on 14/01/20.
//

#ifndef MATHUTILS_INTEGRATION2DTRIANGLE_H
#define MATHUTILS_INTEGRATION2DTRIANGLE_H

#include "Integration2d.h"

namespace mathutils {

  /**
  * Class for computing 2d integrations over a triangle.
  */
  template<typename Tin, typename Tout>
  class Integration2dTriangle : public Integration2d<Tin, Tout> {

   public:

    /// Contructor of the class.
    Integration2dTriangle(Tout (*F)(Tin x), const int &nb_Gauss_points, Vector3d<double> vertex_1,
                          Vector3d<double> vertex_2, Vector3d<double> vertex_3)
        : Integration2d<Tin, Tout>(F, nb_Gauss_points), m_vertex_1(vertex_1), m_vertex_2(vertex_2),
          m_vertex_3(vertex_3) {
    }

    /// This function computes the 2d integration over a triangle.
    virtual Tout Compute() override {

    }

   private:

    /// Three vertices of the triangle.
    Vector3d<double> m_vertex_1;
    Vector3d<double> m_vertex_2;
    Vector3d<double> m_vertex_3;

  };

}

#endif //MATHUTILS_INTEGRATION2DTRIANGLE_H
