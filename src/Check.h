//
// Created by frongere on 16/11/17.
//

#ifndef MATHUTILS_CHECK_H
#define MATHUTILS_CHECK_H

#include <cmath>

namespace mathutils {

    template <class Real=double>
    inline bool IsClose(const Real a, const Real b, const Real rtol = 1e-5, const Real atol = 1e-8) {
        return (fabs(a-b) <= (atol + rtol*fabs(b)));
    }

//    template <class Real=double>
//    inline bool is_unit_vector(chrono::ChVector<Real> vect) {
//        return IsClose(vect.Length2(), 1.);
//    }


}


#endif //MATHUTILS_CHECK_H
