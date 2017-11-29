//
// Created by frongere on 15/11/17.
//

#ifndef FRYDOM_MATRIX_H
#define FRYDOM_MATRIX_H

#include "Eigen/Dense"

#include "Unit.h"

namespace mathutils {


    template <class Scalar=double>
    inline void _GetTrigo(const Scalar pangle, Scalar& cosAngle, Scalar& sinAngle, ANGLE_UNIT unit=RAD) {
        Scalar angle = pangle;
        if (unit == DEG) {
            angle *= MU_PI_180;
        }
        cosAngle = cos(angle);
        sinAngle = sin(angle);
    }


    template <class Scalar=double>
    class Vector2d : public Eigen::Matrix<Scalar, 2, 1> {

    public:

        Vector2d();

        Vector2d(Scalar x, Scalar y);

        void SetNull();

        void Set(Scalar x, Scalar y);

        Vector2d<Scalar> ProjectLocalOnNED(Scalar angle, ANGLE_UNIT unit=RAD) const;
        void ProjectLocalOnNED(Scalar angle, ANGLE_UNIT unit=RAD);

        Vector2d<Scalar> ProjectNEDOnLocal(Scalar angle, ANGLE_UNIT unit=RAD) const;
        void ProjectNEDOnLocal(Scalar angle, ANGLE_UNIT unit=RAD);

        Scalar ProjectOnNorthAxis(Scalar angle, ANGLE_UNIT unit=RAD) const;

        Scalar ProjectOnEastAxis(Scalar angle, ANGLE_UNIT unit=RAD) const;

        Scalar ProjectOnLocalXAxis(Scalar angle, ANGLE_UNIT unit=RAD) const;

        Scalar ProjectOnLocalYAxis(Scalar angle, ANGLE_UNIT unit=RAD) const;

    };

    // =================================================================================================================
    // Functions declarations
    // =================================================================================================================
    template <class Scalar=double>
    Vector2d<Scalar> GetNEDVelocityVector(Scalar plongitudinalSpeed, Scalar plateralSpeed, Scalar angle,
                                          ANGLE_UNIT angleUnit=DEG,
                                          SPEED_UNIT inputSpeedUnit=KNOT, SPEED_UNIT outputSpeedUnit=KNOT);

    template <class Scalar=double>
    Vector2d<Scalar> ProjectLocalOnNED(const Vector2d<Scalar>& vector, Scalar angle, ANGLE_UNIT unit=RAD);

    template <class Scalar=double>
    Vector2d<Scalar> ProjectNEDOnLocal(const Vector2d<Scalar>& vector, Scalar angle, ANGLE_UNIT unit=RAD);

    template <class Scalar=double>
    Scalar ProjectOnNorthAxis(const Vector2d<Scalar>& vector, Scalar angle, ANGLE_UNIT unit=RAD);

    template <class Scalar=double>
    Scalar ProjectOnEastAxis(const Vector2d<Scalar>& vector, Scalar angle, ANGLE_UNIT unit=RAD);

    template <class Scalar=double>
    Scalar ProjectOnLocalXAxis(const Vector2d<Scalar>& vector, Scalar angle, ANGLE_UNIT unit=RAD);

    template <class Scalar=double>
    Scalar ProjectOnLocalYAxis(const Vector2d<Scalar>& vector, Scalar angle, ANGLE_UNIT unit=RAD);

    // =================================================================================================================
    // Vector2d methods implementations
    // =================================================================================================================

    template <class Scalar>
    Vector2d<Scalar>::Vector2d() {
        Eigen::Matrix<Scalar, 2, 1>::setZero();
    }

    template <class Scalar>
    Vector2d<Scalar>::Vector2d(Scalar x, Scalar y) {
        this[0] = x;
        this[1] = y;
    }

    template <class Scalar>
    void Vector2d<Scalar>::Set(Scalar x, Scalar y) {
        this[0] = x;
        this[1] = y;
    }

    template <class Scalar>
    void Vector2d<Scalar>::SetNull() {
        Eigen::Matrix<Scalar, 2, 1>::setZero();
    }

    template <class Scalar>
    Vector2d<Scalar> Vector2d<Scalar>::ProjectLocalOnNED(Scalar angle, ANGLE_UNIT unit) const {
        Scalar c, s;
        _GetTrigo(angle, c, s, unit);
        Vector2d<Scalar> out;

        out.x() = c * this[0] - s * this[1];
        out.y() = s * this[0] + c * this[1];
        return out;
    }

    template <class Scalar>
    void Vector2d<Scalar>::ProjectLocalOnNED(Scalar angle, ANGLE_UNIT unit) {
        Scalar c, s;
        _GetTrigo(angle, c, s, unit);
        this[0] = c * this[0] - s * this[1];
        this[1] = s * this[0] + c * this[1];
    }

    template <class Scalar>
    Vector2d<Scalar> Vector2d<Scalar>::ProjectNEDOnLocal(Scalar angle, ANGLE_UNIT unit) const {
        Scalar c, s;
        _GetTrigo(angle, c, s, unit);
        Vector2d<Scalar> out;
        out.x() =  c * this[0] + s * this[1];
        out.y() = -s * this[0] + c * this[1];
        return out;
    }

    template <class Scalar>
    void Vector2d<Scalar>::ProjectNEDOnLocal(Scalar angle, ANGLE_UNIT unit) {
        Scalar c, s;
        _GetTrigo(angle, c, s, unit);
        this[0] =  c * this[0] + s * this[1];
        this[1] = -s * this[0] + c * this[1];
    }

    template <class Scalar>
    Scalar Vector2d<Scalar>::ProjectOnNorthAxis(Scalar angle, ANGLE_UNIT unit) const {
        Scalar c, s;
        _GetTrigo(angle, c, s, unit);
        return c * this[0] - s * this[1];
    }

    template <class Scalar>
    Scalar Vector2d<Scalar>::ProjectOnEastAxis(Scalar angle, ANGLE_UNIT unit) const {
        Scalar c, s;
        _GetTrigo(angle, c, s, unit);
        return -s * this[0] + c * this[1];
    }

    template <class Scalar>
    Scalar Vector2d<Scalar>::ProjectOnLocalXAxis(Scalar angle, ANGLE_UNIT unit) const {
        Scalar c, s;
        _GetTrigo(angle, c, s, unit);
        return c * this[0] + s * this[1];
    }

    template <class Scalar>
    Scalar Vector2d<Scalar>::ProjectOnLocalYAxis(Scalar angle, ANGLE_UNIT unit) const {
        Scalar c, s;
        _GetTrigo(angle, c, s, unit);
        return -s * this[0] + c * this[1];
    }

    // =================================================================================================================
    // Functions implementations
    // =================================================================================================================

    template <class Scalar>
    Vector2d<Scalar> ProjectLocalOnNED(const Vector2d<Scalar>& vector, Scalar angle, ANGLE_UNIT unit) {
        Vector2d<Scalar> out = vector;
        out.ProjectLocalOnNED(angle, unit);
        return out;
    }

    template <class Scalar>
    Vector2d<Scalar> ProjectNEDOnLocal(const Vector2d<Scalar>& vector, Scalar angle, ANGLE_UNIT unit) {
        Vector2d<Scalar> out = vector;
        out.ProjectNEDOnLocal(angle, unit);
        return out;
    }

    template <class Scalar>
    Scalar ProjectOnNorthAxis(const Vector2d<Scalar>& vector, Scalar angle, ANGLE_UNIT unit) {
        return vector.ProjectOnNorthAxis(angle, unit);
    }

    template <class Scalar>
    Scalar ProjectOnEastAxis(const Vector2d<Scalar>& vector, Scalar angle, ANGLE_UNIT unit) {
        return vector.ProjectOnEastAxis(angle, unit);
    }

    template <class Scalar>
    Scalar ProjectOnLocalXAxis(const Vector2d<Scalar>& vector, Scalar angle, ANGLE_UNIT unit) {
        return vector.ProjectOnLocalXAxis(angle, unit);
    }

    template <class Scalar>
    Scalar ProjectOnLocalYAxis(const Vector2d<Scalar>& vector, Scalar angle, ANGLE_UNIT unit) {
        return vector.ProjectOnLocalYAxis(angle, unit);
    }

    template <class Scalar>
    Vector2d<Scalar> GetNEDVelocityVector(Scalar plongitudinalSpeed, Scalar plateralSpeed, Scalar angle,
                                          ANGLE_UNIT angleUnit,
                                          SPEED_UNIT inputSpeedUnit, SPEED_UNIT outputSpeedUnit) {

        Scalar longitudinalSpeed = plongitudinalSpeed;
        Scalar lateralSpeed      = plateralSpeed;

        if (inputSpeedUnit != outputSpeedUnit) {
            longitudinalSpeed = convert_velocity_unit(longitudinalSpeed, inputSpeedUnit, outputSpeedUnit);
            lateralSpeed      = convert_velocity_unit(lateralSpeed, inputSpeedUnit, outputSpeedUnit);
        }

        auto outVect = Vector2d<Scalar>(longitudinalSpeed, lateralSpeed);
        outVect.ProjectLocalOnNED(angle, angleUnit);  // TODO: verifier qu'on appelle bien la version inline
        return outVect;
    }

}  // end namespace mathutils

#endif //FRYDOM_MATRIX_H
