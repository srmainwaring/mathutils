//
// Created by frongere on 29/12/17.
//

#ifndef MATHUTILS_CONVOLUTION_H
#define MATHUTILS_CONVOLUTION_H

namespace mathutils {


    template <class Real>
    class Convolution {

    private:
        std::vector<Real> m_time;
        std::vector<Real> m_kernel;

    public:
        Convolution() = default;

        explicit Convolution(const std::vector<Real>& time, const std::vector<Real>& kernel)
                : m_time(time), m_kernel(kernel) {}

        void SetKernel(const std::vector<Real>& time, const std::vector<Real>& kernel) {
            m_time = time;
            m_kernel = kernel;
        }



    };



// TODO: creer une classe dans laquelle on peut renseigner le noyau de convolution (reponse impulsionnelle) et fournir un signal
    template <class Real>
    void naiveConvolution(const std::vector<Real>& x, const std::vector<Real>& y) {







        Integrate1d<double> trapz();
    }


}  // end namespace mathutils

#endif //MATHUTILS_CONVOLUTION_H
