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

    template <class Real=double>
    inline bool IsClose(const std::vector<Real>& a, const std::vector<Real>& b, const Real rtol=1e-5, const Real atol=1e-8) {
        assert(a.size() == b.size());
        for (unsigned int i=0; i<a.size(); i++) {
            if (!IsClose(a[i], b[i], rtol, atol)) {
                return false;
            }
        }
        return true;
    }

    template <class Real=double>
    void AssertIsClose(const Real a, const Real b, const Real rtol = 1e-5, const Real atol = 1e-8) {
        assert(IsClose(a, b, rtol, atol));
    }

    template <class Real=double>
    void AssertIsClose(const std::vector<Real>& a, const std::vector<Real>& b, const Real rtol=1e-5, const Real atol=1e-8) {
        assert(IsClose(a, b, rtol, atol));
    }


}


#endif //MATHUTILS_CHECK_H
