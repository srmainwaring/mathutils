//
// Created by frongere on 29/12/17.
//

#ifndef MATHUTILS_CONVOLUTION_H
#define MATHUTILS_CONVOLUTION_H

#include "boost/circular_buffer.hpp"

namespace mathutils {


    template <class Real>
    class Convolution {

    private:
        Real m_timeStep;
        std::vector<Real> m_kernel;

    public:
        Convolution() = default;

        explicit Convolution(const Real dt, const std::vector<Real>& kernel)
                : m_timeStep(dt), m_kernel(kernel) {}

        void SetKernel(const Real& timeStep, const std::vector<Real>& kernel) {
            m_timeStep = timeStep;
            m_kernel = kernel;
        }



    };



    // TODO: creer une classe dans laquelle on peut renseigner le noyau de convolution (reponse impulsionnelle) et fournir un signal
    template <class Real>
    void NaiveConvolution(const std::vector<Real>& kernel, const boost::circular_buffer<Real>& signal, bool reverse=true) {

        unsigned long N = kernel.size();

        assert(signal.size() >= N);

        // Commputing integrand
        std::vector<Real> integrand(N);
        for (unsigned long i=0; i<N; i++) { // TODO: gerer le ccas ou y n'est pas presente sous forme reverse
            integrand[i] = kernel[i] * signal[i];
        }

        // Computing the integral
        Integrate1d<double> trapz();
        // TODO: terminer !!




    }


}  // end namespace mathutils

#endif //MATHUTILS_CONVOLUTION_H
