//
// Created by frongere on 15/11/17.
//

#ifndef FRYDOM_MATRIX_H
#define FRYDOM_MATRIX_H

#include "Eigen/Dense"


namespace mathutils {

    template <class Scalar=double>
    class Vector2d : public Eigen::Matrix<Scalar, 2, 1> {

    public:
        // RotateXY(vector, angle)
        // Vu(vector, angle)
        // Vv(vector, angle)
        // SetBodyVector(pu, pv, pangle)
        // getGivenColumn(matrix, vectorOut&, icol)



    };

}  // end namespace MathUtils

#endif //FRYDOM_MATRIX_H
