//
// Created by frongere on 14/12/17.
//

#ifndef MATHUTILS_MATRIX_H
#define MATHUTILS_MATRIX_H

#include "Eigen/Dense"

namespace mathutils {
    // =================================================================================================================
    // =================================================================================================================
    //                                              DECLARATIONS
    // =================================================================================================================
    // =================================================================================================================


    template <class Scalar>
    class MatrixMN : public Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> {

    public:

        MatrixMN() {};

        MatrixMN(unsigned int nbRows, unsigned int nbCols) : Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>(nbRows, nbCols) {}

        inline Scalar at(unsigned int irow, unsigned int icol) const;

        inline Scalar& at(unsigned int irow, unsigned int icol);


        unsigned int GetNbRows() const;

        unsigned int GetNbCols() const;

        MatrixMN<Scalar> GetDiag() const;

        MatrixMN<Scalar> GetColumn(unsigned int iCol) const;

//        MatrixMN<Scalar>& GetColumn(unsigned int iCol);

        MatrixMN<Scalar> GetRow(unsigned int iRow) const;

//        MatrixMN<Scalar>& GetRow(unsigned int iRow);

        void SetIdentity();

        void SetNull();

        void Transpose();

        void Inverse();

        MatrixMN<Scalar> GetInverse() const;

        MatrixMN<Scalar> GetPseudoInverse(Scalar tol=1e-6) const;

        void Randomize();

        bool IsPositiveSemiDefinite(const Scalar& epsilon=1e-7) const;

        bool IsSymmetric() const;

        bool IsIdentity() const;

        bool IsSquare() const;

        bool IsOrthogonal() const;

        bool IsEqual(const MatrixMN<Scalar>& other, const Scalar& epsilon=1e-12) const;

        void GetQRDecomposition(MatrixMN<Scalar>& Q, MatrixMN<Scalar>& R) const;

        void GetLUDecomposition(MatrixMN<Scalar>& P, MatrixMN<Scalar>& L, MatrixMN<Scalar>& U) const;

        void GetCholeskyDecomposition(MatrixMN<Scalar>& L) const;

        void GetSVDDecomposition(MatrixMN<Scalar>& U, MatrixMN<Scalar>& S, MatrixMN<Scalar>& V) const;




        // ====================================
        // Methods for Eigen inheritance usage
        // ====================================

        // This constructor allows to construct MatrixMN from Eigen expressions
        template <class OtherDerived>
        MatrixMN(const Eigen::MatrixBase<OtherDerived>& other) : Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>(other) {}

        // This method allows to assign Eigen expressions to MatrixMN
        template <class OtherDerived>
        MatrixMN& operator=(const Eigen::MatrixBase<OtherDerived>& other) {
            this->Eigen::Matrix<Scalar, other.rows(), other.cols()>::operator=(other);
            return *this;
        }


    };

    // =================================================================================================================
    // =================================================================================================================
    //                                              FUNCTIONS
    // =================================================================================================================
    // =================================================================================================================

    template <class Scalar>
    MatrixMN<Scalar> Transpose(const MatrixMN<Scalar>& mat) {
        auto tMat = mat;
        tMat.Transpose();
        return tMat;
    }

    template <class Scalar>
    MatrixMN<Scalar> Pinv(const MatrixMN<Scalar>& mat, const Scalar tol=1e-6) {
        return mat.GetPseudoInverse(tol);
    }

    template <class Scalar>
    MatrixMN<Scalar> MakeSymmetricFromUpper(const MatrixMN<Scalar>& other) {
        Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> newMat;
        newMat = other.template selfadjointView<Eigen::Upper>();
        return newMat;
    }

    template <class Scalar>
    MatrixMN<Scalar> MakeSymmetricFromLower(const MatrixMN<Scalar>& other) {
        Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> newMat;
        newMat = other.template selfadjointView<Eigen::Lower>();
        return newMat;
    }





    // =================================================================================================================
    // =================================================================================================================
    //                                              IMPLEMENTATIONS
    // =================================================================================================================
    // =================================================================================================================

    template <class Scalar>
    inline Scalar MatrixMN<Scalar>::at(unsigned int irow, unsigned int icol) const {
        return this->operator()(irow * (this->cols()-1) + icol);
    }

    template <class Scalar>
    inline Scalar& MatrixMN<Scalar>::at(unsigned int irow, unsigned int icol) {
        return this->operator()(irow * (this->cols()-1) + icol);
    }

    template <class Scalar>
    unsigned int MatrixMN<Scalar>::GetNbRows() const {
        return (unsigned int)this->rows();
    }

    template <class Scalar>
    unsigned int MatrixMN<Scalar>::GetNbCols() const {
        return (unsigned int)this->cols();
    }

    template <class Scalar>
    void MatrixMN<Scalar>::SetNull() {
        this->setZero();
    }

    template <class Scalar>
    void MatrixMN<Scalar>::Randomize() {
        this->setRandom();
    }

    template <class Scalar>
    void MatrixMN<Scalar>::SetIdentity() {
        this->setIdentity();
    }

    template <class Scalar>
    MatrixMN<Scalar> MatrixMN<Scalar>::GetColumn(unsigned int iCol) const {
        return this->col(iCol);
    }

    template <class Scalar>
    MatrixMN<Scalar> MatrixMN<Scalar>::GetRow(unsigned int iRow) const {
        return this->row(iRow);
    }

    template <class Scalar>
    MatrixMN<Scalar> MatrixMN<Scalar>::GetDiag() const {
        return this->diagonal();
    }

    template <class Scalar>
    void MatrixMN<Scalar>::Transpose() {
        this->transposeInPlace();
    }

    template <class Scalar>
    void MatrixMN<Scalar>::Inverse() {
        this->swap(this->inverse());
    }

    template <class Scalar>
    MatrixMN<Scalar> MatrixMN<Scalar>::GetInverse() const {
        auto newMat = *this;
        newMat.Inverse();
        return newMat;
    }

    template <class Scalar>
    bool MatrixMN<Scalar>::IsPositiveSemiDefinite(const Scalar& epsilon) const {
        return (this->rows() == this->cols()) && ((this->eigenvalues().real().array() > -epsilon).all());
    }

    template <class Scalar>
    bool MatrixMN<Scalar>::IsSymmetric() const {
        auto a = this->isApprox(this->adjoint());
        // FIXME: finir implementation :::
//        auto a = (this->adjoint() == this);
//        return this
//        if (this->adjoint() == this) {
//            return true;
//        }
        return true;
    }

    template <class Scalar>
    bool MatrixMN<Scalar>::IsIdentity() const {
        return this->isIdentity();
    }

    template <class Scalar>
    bool MatrixMN<Scalar>::IsSquare() const {
        return (this->rows() == this->cols());
    }

    template <class Scalar>
    bool MatrixMN<Scalar>::IsOrthogonal() const {
        return this->isOrthogonal(*this);
    }

    template <class Scalar>
    bool MatrixMN<Scalar>::IsEqual(const MatrixMN<Scalar>& other, const Scalar& epsilon) const {
        return this->isApprox(other, epsilon);
    }

    template <class Scalar>
    void MatrixMN<Scalar>::GetQRDecomposition(MatrixMN<Scalar>& Q, MatrixMN<Scalar>& R) const {
        assert(this->rows() >= this->cols());

        auto QR = this->householderQr();

        // Q
        auto Qfull = QR.householderQ();
        auto thinQ = MatrixMN<Scalar>::Identity(QR.rows(), QR.cols());
        Q = (MatrixMN<Scalar>)(Qfull * thinQ);

        // R
        Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> Rfull;
        Rfull = QR.matrixQR().template triangularView<Eigen::Upper>();
        R = (MatrixMN<Scalar>)(Rfull.block(0, 0, QR.cols(), QR.cols()));  // TODO: verifier que c'est bien cols() cols()...

    }

    template <class Scalar>
    void MatrixMN<Scalar>::GetLUDecomposition(MatrixMN<Scalar>& P, MatrixMN<Scalar>& L, MatrixMN<Scalar>& U) const {

        // PA = LU

        assert(this->rows() == this->cols());

        auto LU = this->partialPivLu();

        // Permutation matrix
        Eigen::Matrix<unsigned int, Eigen::Dynamic, Eigen::Dynamic> PP;
        PP = LU.permutationP();
        Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> PPScalar = PP.template cast<Scalar>();
        P = (MatrixMN<Scalar>)(PPScalar);

        // L
        Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> Lpartial;
        Lpartial = LU.matrixLU().template triangularView<Eigen::StrictlyLower>();
        L = (MatrixMN<Scalar>)(Lpartial + Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>::Identity(LU.rows(), LU.cols()));

        // U
        Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> Utmp;
        Utmp = LU.matrixLU().template triangularView<Eigen::Upper>();
        U = (MatrixMN<Scalar>)(Utmp);

    }

    template <class Scalar>
    void MatrixMN<Scalar>::GetCholeskyDecomposition(MatrixMN<Scalar>& L) const {

        // Verifying that the matrix is square
        if (!IsSquare()) {
            throw std::runtime_error("In Cholesky decomposition, matrix must be square");
        }

        // Verifying that the matrix is symmetric
        if (!IsSymmetric()) {
            throw std::runtime_error("In Cholesky decomposition, matrix must be self adjoint (symmetric if real coefficients)");
        }

        auto CHOL = this->llt();

        if (CHOL.info() == Eigen::NumericalIssue) {
            throw std::runtime_error("In Cholesky decomposition, matrix must be positive semi definite and this one is possibly not");
        }

        Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> Ltmp;
        Ltmp = CHOL.matrixL();
        L = (MatrixMN<Scalar>)(Ltmp);
    }

    template <class Scalar>
    void MatrixMN<Scalar>::GetSVDDecomposition(MatrixMN<Scalar>& U, MatrixMN<Scalar>& S, MatrixMN<Scalar>& V) const {

        auto SVD = this->jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV);

        Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> Utmp, Stmp, Vtmp;

        // U
        Utmp = SVD.matrixU();
        U = (MatrixMN<Scalar>)(Utmp);

        // S
        Stmp = SVD.singularValues();
        S = (MatrixMN<Scalar>)(Stmp);

        // V
        Vtmp = SVD.matrixV();
        V = (MatrixMN<Scalar>)(Vtmp);
    }

    template <class Scalar>
    MatrixMN<Scalar> MatrixMN<Scalar>::GetPseudoInverse(const Scalar tol) const {

        // SVD decomposition
        MatrixMN<Scalar> U, S, V;
        GetSVDDecomposition(U, S, V);

        for (unsigned int i=0; i < S.rows(); ++i) {
            if (S(i) > tol) {
                S(i) = 1. / S(i);
            } else {
                S(i) = 0.;
            }
        }

        return V * S.asDiagonal() * U.adjoint();
    }


}  // end namespace mathutils

#endif //MATHUTILS_MATRIX_H
