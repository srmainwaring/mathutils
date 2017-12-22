//
// Created by frongere on 22/12/17.
//

#ifndef MATHUTILS_FFT_H
#define MATHUTILS_FFT_H

#include "kiss_fft.h"


namespace mathutils {

    template <class Real>
    class FFT {

    private:
        kiss_fft_cfg     m_fft_cfg;

//        unsigned int NFFT;
//        bool Inverse;


    public:
        FFT() {

//            kiss_fft_alloc();

        }

        FFT(unsigned int NFFT) {
            SetNbPoints(NFFT);
        }

        ~FFT() {
            std::cout << "Destroying FFT object" << std::endl;
            kiss_fft_free(m_fft_cfg);
        }

        void SetNbPoints(unsigned int NFFT) {
            m_fft_cfg = kiss_fft_alloc(NFFT, 0, 0, 0);
        }



        std::vector<Real> fft(const std::vector<Real>& f) {

        }





    };




}  // end namespace mathutils


#endif //MATHUTILS_FFT_H
