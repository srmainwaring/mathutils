//
// Created by frongere on 15/11/17.
//

#ifndef FRYDOM_MATRIX_H
#define FRYDOM_MATRIX_H

#include "Eigen/Dense"

#include "Unit.h"
#include "Angles.h"

namespace mathutils {
    // See the following link for Eigen::Matrix inheritannce :
    // eigen.tuxfamily.org/dox/TopicCustomizing_InheritingMatrix.html




    template <class Scalar=double>
    class Vector2d : public Eigen::Matrix<Scalar, 2, 1> {

    public:

        Vector2d();

        Vector2d(Scalar x, Scalar y);

        // This constructor allows to construct Vector2d from Eigen expressions
        template <class OtherDerived>
        Vector2d(const Eigen::MatrixBase<OtherDerived>& other) : Eigen::Matrix<Scalar, 2, 1>(other) {}

        // This method allows to assign Eigen expressions to Vector2d
        template <class OtherDerived>
        Vector2d& operator=(const Eigen::MatrixBase<OtherDerived>& other) {
            this->Eigen::Matrix<Scalar, 2, 1>::operator=(other);
            return *this;
        }

        inline Scalar at(unsigned int index) const { return this->operator[](index); }

        inline Scalar& at(unsigned int index) { return this->operator[](index); }

        inline void SetNull();

        inline void Set(Scalar x, Scalar y);

        inline void ProjectLocalOnNED(Vector2d<Scalar>& out, Scalar angle, ANGLE_UNIT unit=RAD) const;
        inline void ProjectLocalOnNED(Scalar angle, ANGLE_UNIT unit=RAD);  // TODO: voir pour le inplace

        inline void ProjectNEDOnLocal(Vector2d<Scalar>& out, Scalar angle, ANGLE_UNIT unit=RAD) const;
        inline void ProjectNEDOnLocal(Scalar angle, ANGLE_UNIT unit=RAD);  // TODO: voir pour le inplace

        inline Scalar ProjectOnNorthAxis(Scalar angle, ANGLE_UNIT unit=RAD) const;

        inline Scalar ProjectOnEastAxis(Scalar angle, ANGLE_UNIT unit=RAD) const;

        inline Scalar ProjectOnLocalXAxis(Scalar angle, ANGLE_UNIT unit=RAD) const;

        inline Scalar ProjectOnLocalYAxis(Scalar angle, ANGLE_UNIT unit=RAD) const;

        inline Scalar infNorm() const;

    };

    // =================================================================================================================
    // Functions declarations
    // =================================================================================================================
    template <class Scalar=double>
    Vector2d<Scalar> GetNEDVelocityVector(Scalar plongitudinalSpeed, Scalar plateralSpeed, Scalar angle,
                                          ANGLE_UNIT angleUnit=RAD,
                                          SPEED_UNIT inputSpeedUnit=KNOT, SPEED_UNIT outputSpeedUnit=MS);

    template <class Scalar>
    Vector2d<Scalar> ProjectLocalOnNED(const Vector2d<Scalar>& vector, Scalar angle, ANGLE_UNIT unit=RAD);

    template <class Scalar>
    Vector2d<Scalar> ProjectNEDOnLocal(const Vector2d<Scalar>& vector, Scalar angle, ANGLE_UNIT unit=RAD);

    template <class Scalar>
    Scalar ProjectOnNorthAxis(const Vector2d<Scalar>& vector, Scalar angle, ANGLE_UNIT unit=RAD);

    template <class Scalar>
    Scalar ProjectOnEastAxis(const Vector2d<Scalar>& vector, Scalar angle, ANGLE_UNIT unit=RAD);

    template <class Scalar>
    Scalar ProjectOnLocalXAxis(const Vector2d<Scalar>& vector, Scalar angle, ANGLE_UNIT unit=RAD);

    template <class Scalar>
    Scalar ProjectOnLocalYAxis(const Vector2d<Scalar>& vector, Scalar angle, ANGLE_UNIT unit=RAD);

    // =================================================================================================================
    // Vector2d methods implementations
    // =================================================================================================================

    template <class Scalar>
    Vector2d<Scalar>::Vector2d() : Eigen::Matrix<Scalar, 2, 1>() {
        Eigen::Matrix<Scalar, 2, 1>::setZero();
    }

    template <class Scalar>
    Vector2d<Scalar>::Vector2d(Scalar x, Scalar y) {
        this->operator[](0) = x;
        this->operator[](1) = y;
    }

    template <class Scalar>
    void Vector2d<Scalar>::Set(Scalar x, Scalar y) {
        this->operator[](0) = x;
        this->operator[](1) = y;
    }

    template <class Scalar>
    void Vector2d<Scalar>::SetNull() {
        Eigen::Matrix<Scalar, 2, 1>::setZero();
    }

    template <class Scalar>
    void Vector2d<Scalar>::ProjectLocalOnNED(Vector2d<Scalar>& out, Scalar angle, ANGLE_UNIT unit) const {
        Scalar c, s;
        _GetTrigo(angle, c, s, unit);

        out[0] = c * at(0) - s * at(1);
        out[1] = s * at(0) + c * at(1);  // TODO: utiliser l'autre methode
    }

    template <class Scalar>
    void Vector2d<Scalar>::ProjectLocalOnNED(Scalar angle, ANGLE_UNIT unit) {
        Scalar c, s;
        _GetTrigo(angle, c, s, unit);
        auto x = c * at(0) - s * at(1);
        at(1) = s * at(0) + c * at(1);
        at(0) = x;
    }

    template <class Scalar>
    void Vector2d<Scalar>::ProjectNEDOnLocal(Vector2d<Scalar>& out, Scalar angle, ANGLE_UNIT unit) const {
        Scalar c, s;
        _GetTrigo(angle, c, s, unit);
        out.x() =  c * at(0) + s * at(1);
        out.y() = -s * at(0) + c * at(1);  // TODO: utiliser l'autre methode
    }

    template <class Scalar>
    void Vector2d<Scalar>::ProjectNEDOnLocal(Scalar angle, ANGLE_UNIT unit) {
        Scalar c, s;
        _GetTrigo(angle, c, s, unit);
        auto x =  c * at(0) + s * at(1);
        at(1) = -s * at(0) + c * at(1);
        at(0) = x;
    }

    template <class Scalar>
    Scalar Vector2d<Scalar>::ProjectOnNorthAxis(Scalar angle, ANGLE_UNIT unit) const {
        Scalar c, s;
        _GetTrigo(angle, c, s, unit);
        return c * at(0) - s * at(1);
    }

    template <class Scalar>
    Scalar Vector2d<Scalar>::ProjectOnEastAxis(Scalar angle, ANGLE_UNIT unit) const {
        Scalar c, s;
        _GetTrigo(angle, c, s, unit);
        return s * at(0) + c * at(1);
    }

    template <class Scalar>
    Scalar Vector2d<Scalar>::ProjectOnLocalXAxis(Scalar angle, ANGLE_UNIT unit) const {
        Scalar c, s;
        _GetTrigo(angle, c, s, unit);
        return c * at(0) + s * at(1);
    }

    template <class Scalar>
    Scalar Vector2d<Scalar>::ProjectOnLocalYAxis(Scalar angle, ANGLE_UNIT unit) const {
        Scalar c, s;
        _GetTrigo(angle, c, s, unit);
        return -s * at(0) + c * at(1);
    }

    template <class Scalar>
    Scalar Vector2d<Scalar>::infNorm() const {
        return this->Eigen::Matrix<Scalar, 2, 1>::maxCoeff();;
    }

    // =================================================================================================================
    // Functions implementations
    // =================================================================================================================

    template<class Scalar>
    Vector2d<Scalar> ProjectLocalOnNED(const Vector2d<Scalar> &vector, Scalar angle, ANGLE_UNIT unit) {
        Vector2d<Scalar> out = vector;
        out.ProjectLocalOnNED(angle, unit);
        return out;
    }

    template<class Scalar>
    Scalar ProjectOnNorthAxis(const Vector2d<Scalar> &vector, Scalar angle, ANGLE_UNIT unit) {
        return vector.ProjectOnNorthAxis(angle, unit);
    }

    template<class Scalar>
    Vector2d<Scalar> ProjectNEDOnLocal(const Vector2d<Scalar> &vector, Scalar angle, ANGLE_UNIT unit) {
        Vector2d<Scalar> out = vector;
        out.ProjectNEDOnLocal(angle, unit);
        return out;
    }

    template<class Scalar>
    Scalar ProjectOnEastAxis(const Vector2d<Scalar> &vector, Scalar angle, ANGLE_UNIT unit) {
        return vector.ProjectOnEastAxis(angle, unit);
    }

    template<class Scalar>
    Scalar ProjectOnLocalXAxis(const Vector2d<Scalar> &vector, Scalar angle, ANGLE_UNIT unit) {
        return vector.ProjectOnLocalXAxis(angle, unit);
    }

    template<class Scalar>
    Scalar ProjectOnLocalYAxis(const Vector2d<Scalar> &vector, Scalar angle, ANGLE_UNIT unit) {
        return vector.ProjectOnLocalYAxis(angle, unit);
    }

    template<class Scalar=double>
    Vector2d<Scalar>
    GetNEDVelocityVector(Scalar plongitudinalSpeed, Scalar plateralSpeed, Scalar angle, ANGLE_UNIT angleUnit,
                         SPEED_UNIT inputSpeedUnit, SPEED_UNIT outputSpeedUnit) {

        Scalar longitudinalSpeed = plongitudinalSpeed;
        Scalar lateralSpeed      = plateralSpeed;

        if (inputSpeedUnit != outputSpeedUnit) {
            longitudinalSpeed = convert_velocity_unit(longitudinalSpeed, inputSpeedUnit, outputSpeedUnit);
            lateralSpeed      = convert_velocity_unit(lateralSpeed, inputSpeedUnit, outputSpeedUnit);
        }

        auto localVect = Vector2d<Scalar>(longitudinalSpeed, lateralSpeed);
        return ProjectLocalOnNED(localVect, angle, angleUnit);
//        return outVect.ProjectLocalOnNED(angle, angleUnit);  // TODO: verifier qu'on appelle bien la version inline
    }

}  // end namespace mathutils

#endif //FRYDOM_MATRIX_H
