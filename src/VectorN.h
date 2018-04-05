//
// Created by frongere on 25/01/18.
//

#ifndef MATHUTILS_VECTORN_H
#define MATHUTILS_VECTORN_H

#include "Eigen/Dense"
#include "Matrix.h"

namespace mathutils {

    // TODO: travailler sur les operateurs + et *
    // TODO: Rajouter methode Resize()


    template <class Scalar>
    class VectorN : public Eigen::Matrix<Scalar, Eigen::Dynamic, 1> {

    public:

        VectorN() {}

        VectorN(unsigned int nbRows) : Eigen::Matrix<Scalar, Eigen::Dynamic, 1>(nbRows) {}


        inline Scalar at(unsigned int n) const;

        inline Scalar& at(unsigned int n);

        void Randomize();

        void SetNull();

        void Sort(bool ascending=true);

        void Reverse();

        bool IsEqual(const VectorN<Scalar>& other, const Scalar& epsilon=1e-12);


        MatrixMN<Scalar> GetMatrixSquare() const;



        // =====================================================================
        // Methods for Eigen inheritance usage
        // =====================================================================

        // This constructor allows to construct VectorN from Eigen expressions
        template <class OtherDerived>
        VectorN(const Eigen::MatrixBase<OtherDerived>& other) : Eigen::Matrix<Scalar, Eigen::Dynamic, 1>(other) {}

        // This method allows to assign Eigen expressions to VectorN
        template <class OtherDerived>
        VectorN& operator=(const Eigen::MatrixBase<OtherDerived>& other) {
            this->Eigen::Matrix<Scalar, other.rows(), 1>::operator=(other);
            return *this;
        }



    };


    // =================================================================================================================
    // =================================================================================================================
    //                                              IMPLEMENTATIONS
    // =================================================================================================================
    // =================================================================================================================

    template <class Scalar>
    inline Scalar VectorN<Scalar>::at(unsigned int n) const {
        return this->operator()(n);
    }

    template <class Scalar>
    inline Scalar& VectorN<Scalar>::at(unsigned int n) {
        return this->operator()(n);
    }

    template <class Scalar>
    void VectorN<Scalar>::Randomize() {
        this->setRandom();
    }

    template <class Scalar>
    void VectorN<Scalar>::SetNull() {
        this->setZero();
    }

    template <class Scalar>
    void VectorN<Scalar>::Sort(bool ascending) {
        std::sort(this->data(), this->data()+this->size());

        if (!ascending) Reverse();

    }

    MatrixMN<Scalar> VectorN<Scalar>::GetMatrixSquare() const
    {
        MatrixMN<Scalar> Res(this->size(),this->size());
        VectorN<Scalar> V;
        Res = this * this->transpose();
        return Res;
    }

    template <class Scalar>
    void VectorN<Scalar>::Reverse() {
        std::reverse(this->data(), this->data()+this->size());
    }

    template <class Scalar>
    bool VectorN<Scalar>::IsEqual(const VectorN<Scalar>& other, const Scalar& epsilon) {
        return this->isApprox(other, epsilon);
    }

}  // End namespace mathutils



#endif //MATHUTILS_VECTORN_H
