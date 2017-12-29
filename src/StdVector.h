//
// Created by frongere on 29/12/17.
//

#ifndef MATHUTILS_STDVECTOR_H
#define MATHUTILS_STDVECTOR_H

namespace mathutils {

    template <class Scalar>
    inline std::vector<Scalar> Add(const std::vector<Scalar>& v1, const std::vector<Scalar>& v2) {
        unsigned int n = v1.size();
        assert(v2.size() == n);

        std::vector<Scalar> out(n);
        for (unsigned int i=0; i<n; i++) {
            out[i] = v1[i] + v2[i];
        }
        return out;
    }

    template <class Scalar>
    inline std::vector<Scalar> Minus(const std::vector<Scalar>& v1, const std::vector<Scalar>& v2) {
        unsigned int n = v1.size();
        assert(v2.size() == n);

        std::vector<Scalar> out(n);
        for (unsigned int i=0; i<n; i++) {
            out[i] = v1[i] - v2[i];
        }
        return out;
    }

    template <class Scalar>
    inline std::vector<Scalar> Mult(const std::vector<Scalar>& v1, const std::vector<Scalar>& v2) {
        unsigned int n = v1.size();
        assert(v2.size() == n);

        std::vector<Scalar> out(n);
        for (unsigned int i=0; i<n; i++) {
            out[i] = v1[i] * v2[i];
        }
        return out;
    }

    template <class Scalar>
    inline void Scale(std::vector<Scalar>& v, const Scalar alpha) {
        for (Scalar& vi: v) {
            vi *= alpha;
        }
    }

    template <class Scalar>
    inline void Pow(std::vector<Scalar>& v, const Scalar p) {
        for (Scalar& vi: v) {
            vi = pow(vi, p);
        }
    }

    template <class Scalar>
    inline void Negate(std::vector<Scalar>& v) {
        Scale(v, -1.);
    }

}


#endif //MATHUTILS_STDVECTOR_H
