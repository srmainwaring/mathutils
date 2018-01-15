//
// Created by frongere on 29/11/17.
//

#ifndef MATHUTILS_MATHS_H
#define MATHUTILS_MATHS_H

namespace mathutils {

    template <typename T>
    int sgn(T val) {
        return (T(0) < val) - (val < T(0));
    }

}  // end namespace mathutils


#endif //MATHUTILS_MATHS_H
