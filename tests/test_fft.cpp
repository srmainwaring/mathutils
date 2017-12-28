//
// Created by frongere on 22/12/17.
//

#include "MathUtils.h"

using namespace mathutils;


// TODO: a mettre dans FFT.h
unsigned int NextPow2(const unsigned int n, const bool returnValue=true) {
    unsigned int p = 0;
    unsigned int val;
    while (true) {
        val = (unsigned int)pow(2, p);
        if (val >= n) {
            if (returnValue) {
                return val;
            } else {
                return p;
            }
        }
        p++;
    }
}




template <class Scalar>
std::vector<Scalar> MakeDataSin(const std::vector<Scalar>& timeVector) {

    auto len = timeVector.size();

    std::vector<Scalar> vector;
    vector.reserve(len);

    for (auto t: timeVector) {
        vector.push_back(sin(t));
    }

    return vector;

}


int main(int argc, char* argv[]) {


    // Building data
    auto time = linspace<double>(0, 2.*MU_2PI, NextPow2(1000));

    auto timeVector = MakeDataSin<double>(time);


    // Building the fft object
    FFT<double> fft;
    fft.ScalingON();
    fft.HalfSpectrumON();

    // Computing FFT
    std::vector<std::complex<double>> freqVect;
    fft.fft(freqVect, timeVector);

    // Computing IFFT from spectrum
    std::vector<double> timeVectRebuild;
    fft.ifft(timeVectRebuild, freqVect);

    AssertIsClose(timeVector, timeVectRebuild);

    return 0;
}