//
// Created by frongere on 22/12/17.
//

#ifndef MATHUTILS_FFT_H
#define MATHUTILS_FFT_H

//#include "kiss_fft.h"

#include <unsupported/Eigen/FFT>

namespace mathutils {

    template <class Scalar>
    class FFT : protected Eigen::FFT<Scalar> {

    public:

        void fft(std::vector<std::complex<Scalar>>& freqVector, const std::vector<Scalar>& vector) {
            this->fwd(freqVector, vector);
        }

        void ifft(std::vector<Scalar>& vector, const std::vector<std::complex<Scalar>>& freqVector) {
            this->inv(vector, freqVector);
        }

        void HalfSpectrumON() {
            this->SetFlag(Eigen::FFT<Scalar>::HalfSpectrum);
//            this->HalfSpectrum = true;
        }

        void HalfSpectrumOFF() {
            this->ClearFlag(Eigen::FFT<Scalar>::HalfSpectrum);
//            this->HalfSpectrum = false;
        }

        void ScalingON() {
            this->ClearFlag(Eigen::FFT<Scalar>::Unscaled);
//            this->Unscaled = false;
        }

        void ScalingOFF() {
            this->SetFlag(Eigen::FFT<Scalar>::Unscaled);
//            this->Unscaled = true;
        }

    };








//    template <class Real>
//    class FFT {
//
//    private:
//        kiss_fft_cfg     m_fft_cfg;
//
////        unsigned int NFFT;
////        bool Inverse;
//
//
//    public:
//        FFT() {
//
////            kiss_fft_alloc();
//
//        }
//
//        FFT(unsigned int NFFT) {
//            SetNbPoints(NFFT);
//        }
//
//        ~FFT() {
//            std::cout << "Destroying FFT object" << std::endl;
//            kiss_fft_free(m_fft_cfg);
//        }
//
//        void SetNbPoints(unsigned int NFFT) {
//            m_fft_cfg = kiss_fft_alloc(NFFT, 0, 0, 0);
//        }
//
//
//
//        std::vector<Real> fft(const std::vector<Real>& f) {
//
//        }
//
//
//
//
//
//    };




}  // end namespace mathutils


#endif //MATHUTILS_FFT_H
