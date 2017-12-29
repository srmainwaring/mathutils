//
// Created by frongere on 28/12/17.
//

#ifndef MATHUTILS_SIGNALS_H
#define MATHUTILS_SIGNALS_H

#include "StdVector.h"

namespace mathutils {

    template <class Real>
    std::vector<Real> Sinc(const std::vector<Real>& x) {
        unsigned long n = x.size();
        std::vector<Real> out(n);
        for (unsigned long i=0; i<n; i++) {
            if (x[i] == 0.) {
                out[i] = 1.;
            } else {
                out[i] = sin(x[i]) / x[i];
            }
        }
        return out;
    }

    template <class Real>
    std::vector<Real> Sinc(const Real xmin, const Real xmax, const unsigned int n) {
        auto x = linspace<Real>(xmin, xmax, n);
        return Sinc(x);
    }

    template <class Real>
    std::vector<Real> Heaviside(const std::vector<Real>& x, const Real x_step) {
        unsigned long n = x.size();

        std::vector<Real> out(n);
        Real xval;
        for (unsigned long i=0; i<n; i++) {
            xval = x[i];
            if (xval < x_step) {
                out[i] = 0.;
            } else if (xval > x_step) {
                out[i] = 1.;
            } else {
                out[i] = 0.5;
            }
        }
        return out;
    }

    template <class Real>
    std::vector<Real> Door(const std::vector<Real>& x, const Real x_min, const Real x_max, const Real height=1.) {
        auto out = Heaviside(x, x_min);
        out = Minus(out, Heaviside<Real>(x, x_max));
        Scale<Real>(out, height);
        return out;
    }

    template <class Real>
    std::vector<Real> Dirac(const std::vector<Real>& x, const Real x_impulse, const Real height=1.) {
        unsigned long n = x.size();

        Real dx_2 = 0.5 * (x[1] - x[0]);

        std::vector<Real> out(n);
        for (unsigned long i=0; i<n; i++) {
            if (abs(x[i] - x_impulse) < dx_2) { // FIXME: le critere n'est pas robuste
                out[i] = height;
            } else {
                out[i] = 0.;
            }
        }
        return out;
    }


//    template <class Real>
//    std::vector<Real> Heaviside(const Real x, const Real xmin, const Real xmax, const unsigned int n) {
//        assert(x >= xmin && w <= xmax);
//
//        auto xVect = linspace<Real>(xmin, xmax, n);
//
//        auto out = std::vector<Real>(n);
//
//        for (unsigned int i=0; i<n; i++) {
//            if (xVect[i] < x) {
//                out[i] = 0.;
//            } else {
//
//            }
//
//        }
//        return out;
//    }


}  // end namespace mathutils


#endif //MATHUTILS_SIGNALS_H
