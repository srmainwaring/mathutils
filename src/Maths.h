//
// Created by frongere on 29/11/17.
//

#ifndef FRYDOM_MATHS_H
#define FRYDOM_MATHS_H

namespace mathutils {

    template <typename T>
    int sgn(T val) {
        return (T(0) < val) - (val < T(0));
    }

}  // end namespace mathutils


#endif //FRYDOM_MATHS_H
