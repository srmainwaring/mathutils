//
// Created by frongere on 29/11/17.
//

#ifndef MATHUTILS_VECTOR3D_H
#define MATHUTILS_VECTOR3D_H

#include "Eigen/Dense"

#include "Unit.h"
#include "Matrix.h"

namespace mathutils {
    // See the following link for Eigen::Matrix inheritance :
    // eigen.tuxfamily.org/dox/TopicCustomizing_InheritingMatrix.html

    // =================================================================================================================
    // =================================================================================================================
    //                                              DECLARATIONS
    // =================================================================================================================
    // =================================================================================================================

    template <class Scalar=double>
    class Vector3d : public Eigen::Matrix<Scalar, 3, 1> {

    public:

        Vector3d();

        Vector3d(Scalar x, Scalar y, Scalar z);

        // This constructor allows to construct Vector2d from Eigen expressions
        template <class OtherDerived>
        Vector3d(const Eigen::MatrixBase<OtherDerived>& other) : Eigen::Matrix<Scalar, 3, 1>(other) {}

        // This method allows to assign Eigen expressions to Vector3d
        template <class OtherDerived>
        Vector3d& operator=(const Eigen::MatrixBase<OtherDerived>& other) {
            this->Eigen::Matrix<Scalar, 3, 1>::operator=(other);
            return *this;
        }

        inline Scalar at(unsigned int index) const { return this->operator[](index); }

        inline Scalar& at(unsigned int index) { return this->operator[](index); }

        inline void SetNull();

        inline void Set(Scalar x, Scalar y, Scalar z);

        inline Scalar infNorm() const;

        void print() const;

    };

    // =================================================================================================================
    // Functions declarations
    // =================================================================================================================





    // =================================================================================================================
    // =================================================================================================================
    //                                              IMPLEMENTATIONS
    // =================================================================================================================
    // =================================================================================================================

    // =================================================================================================================
    // Vector3d methods implementations
    // =================================================================================================================

    template <class Scalar>
    Vector3d<Scalar>::Vector3d() : Eigen::Matrix<Scalar, 3, 1>() {
        SetNull();
    }

    template <class Scalar>
    Vector3d<Scalar>::Vector3d(Scalar x, Scalar y, Scalar z) {
        Set(x, y, z);
    }

    template <class Scalar>
    void Vector3d<Scalar>::SetNull() {
        Eigen::Matrix<Scalar, 3, 1>::setZero();
    }

    template <class Scalar>
    void Vector3d<Scalar>::Set(Scalar x, Scalar y, Scalar z) {
        this->operator[](0) = x;
        this->operator[](1) = y;
        this->operator[](2) = z;
    }

    template <class Scalar>
    Scalar Vector3d<Scalar>::infNorm() const {
        return this->Eigen::Matrix<Scalar, 3, 1>::template lpNorm<Eigen::Infinity>();
    }

    template <class Scalar>
    void Vector3d<Scalar>::print() const {
        std::cout << std::endl << *this << std::endl;
    }

    // =================================================================================================================
    // Functions implementations
    // =================================================================================================================

    template <class Scalar>
    MatrixMN<Scalar> outer(const Vector3d<Scalar>& v1, const Vector3d<Scalar>& v2) {
        Scalar a, b, c;
        a = v1(0);
        b = v1(1);
        c = v1(2);

        Scalar x, y, z;
        x = v2(0);
        y = v2(1);
        z = v2(2);

        MatrixMN<Scalar> mat(3, 3);
        mat.at(0, 0) = a * x;
        mat.at(0, 1) = a * y;
        mat.at(0, 2) = a * z;
        mat.at(1, 0) = b * x;
        mat.at(1, 1) = b * y;
        mat.at(1, 2) = b * z;
        mat.at(2, 0) = c * x;
        mat.at(2, 1) = c * y;
        mat.at(2, 2) = c * z;

        return mat;
    }

    template <class Scalar>
    MatrixMN<Scalar> outer(const Vector3d<Scalar>& v) {
        return outer(v, v);
    }

    template <class Scalar>
    double inner(const Vector3d<Scalar>& v1, const Vector3d<Scalar>& v2) {
        return v1(0)*v2(0) + v1(1)*v2(1) + v1(2)*v2(2);
    }

}  // end namespace mathutils


#endif //MATHUTILS_VECTOR3D_H
