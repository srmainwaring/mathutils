//
// Created by frongere on 16/11/17.
//

#ifndef MATHUTILS_UNIT_H
#define MATHUTILS_UNIT_H

// TODO: voir à intégrer le framework de gestion des unités :   github.com/nholthaus/units

#include <iostream>  // TODO: retirer et mettre en place un systeme d'erreur partout ou c'est utilise...
#include "Constants.h"


/*
 * Here we define unit related functions such as conversions between different units
 */

namespace mathutils {


    /* ===========================================================================================
     * ANGLES Conversions functions
     */

    enum ANGLE_UNIT {
        DEG,
        RAD
    };

    /// CONVERSION DEG->RAD
    template <class T>
    inline T radians(const T a) {
        return a * MU_PI_180;
    }

    /// CONVERSION RAD->DEG
    template <class T>
    inline T degrees(const T a) {
        return a / MU_PI_180;
    }

    // TODO: both two following functions should be realized as MACROS for performance
    /// CONVERSION DEG->RAD
    template <class Real=double>
    inline Real deg2rad(const Real a) {
        return a * MU_PI_180;
    }

    /// CONVERSION RAD->DEG
    template <class Real=double>
    inline Real rad2deg(const Real a) {
        return a / MU_PI_180;
    }

    /* ===========================================================================================
     * Velocities conversions functions
     */

    enum SPEED_UNIT {
        MS,     ///> M/S
        KNOT,   ///> NAUTICAL KNOTS
        KMH     ///> KM/H
    };


    /// Convert nautical knots into m/s
    template <class Real=double>
    inline Real KNOT2MS(Real velocity) {
        return velocity * MU_KNOT;
    }

    /// Convert m/s into nautical knots
    template <class Real=double>
    inline Real MS2KNOT(Real velocity) {
        return velocity / MU_KNOT;
    }

    /// Convert km/h into m/s
    template <class Real=double>
    inline Real KMH2MS(Real velocity) {
        return velocity * MU_KMH;
    }

    /// Convert m/s into km/h
    template <class Real=double>
    inline Real MS2KMH(Real velocity) {
        return velocity / MU_KMH;
    }

    /// Convert km/h into nautical knots
    template <class Real=double>
    inline Real KMH2KNOT(Real velocity) {
        return MS2KNOT(KMH2MS(velocity));
    }

    /// Convert nautical knots into km/h
    template <class Real=double>
    inline Real KNOT2KMH(Real velocity) {
        return MS2KMH(KNOT2MS(velocity));
    }

    /// Converts a velocity from a unit to another
    template <class Real=double>
    inline Real convert_velocity_unit(const Real velocity, SPEED_UNIT current_unit, SPEED_UNIT new_unit= MS) {

        Real new_vel;
        // Expressing in M/S
        switch (current_unit) {
            case MS:
                new_vel = velocity;
                break;
            case KMH:
                new_vel = KMH2MS(velocity);
                break;
            case KNOT:
                new_vel = KNOT2MS(velocity);
        }

        // EXPRESSING IN NEW UNIT
        switch (new_unit) {
            case MS:
                return new_vel;
            case KMH:
                return MS2KMH(new_vel);
            case KNOT:
                return MS2KNOT(new_vel);
        }

    }


    /* ===========================================================================================
     * Frequency conversion functions
     */

    enum FREQUENCY_UNIT {
        HZ,     ///< Hertz (s**-1)
        RADS,   ///< rad/s
        S,      ///< seconds (period)
        RPM     ///< Round Per Minute
    };  // TODO: ajouter la conversion avec RPM ...

    template <class Real=double>
    inline Real HZ2RADS(const Real hz) {
        return MU_2PI * hz;
    }

    template <class Real=double>
    inline Real RADS2HZ(const Real rads) {
        return rads / MU_2PI;
    }

    template <class Real=double>
    inline Real HZ2S(const Real hz) {
        return 1. / hz;
    }

    template <class Real=double>
    inline Real S2HZ(const Real s) {
        return 1. / s;
    }

    template <class Real=double>
    inline Real RADS2S(const Real rads) {
        return MU_2PI / rads;
    }

    template <class Real=double>
    inline Real S2RADS(const Real s) {
        return MU_2PI / s;
    }

    template <class Real=double>
    inline Real convert_frequency(const Real in, FREQUENCY_UNIT src_unit, FREQUENCY_UNIT target_unit) {

        if (src_unit == target_unit) return in;

        Real piv;  // Must be in Hz

        switch (src_unit) {
            case HZ:
                piv = in;
                break;
            case RADS:
                piv = RADS2HZ(in);
                break;
            case S:
                piv = S2HZ(in);
                break;
        }

        switch (target_unit) {
            case RADS:
                return HZ2RADS(piv);
            case S:
                return HZ2S(piv);
            case HZ:
                std::cout << "Impossible case in frequency conversion !!" << std::endl;
                break;
        }

    }


}

#endif //MATHUTILS_UNIT_H
