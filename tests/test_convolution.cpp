//
// Created by frongere on 29/12/17.
//

#include "MathUtils.h"
#include "matplotlibcpp.h"


using namespace mathutils;

template <class Real>
void ZeroPad(std::vector<Real>& vector, unsigned int n) {
    assert(n >= vector.size());

    vector.reserve(n);
    for (unsigned int i=0; i<vector.size()-n; i++) {
        vector.push_back(0.);
    }
}




int main(int argc, char* argv[]) {

    unsigned int n = 500;
    unsigned int m = 700;

    // Creating the impulse response signal
    auto tr = linspace<double>(0, 100, n);

    std::vector<double> r(n);
    r[0] = 1.;
    for (uint i=1; i<n; i++) {
        r[i] = sin(tr[i]) / tr[i];
    }

    // Creating the convolution engine
    Convolution<double> conv(tr, r);


    // Creating the signal
    auto ts = linspace<double>(0., (double)m*100/(double)n, m);
    std::vector<double> s(m);
    for (uint i=0; i<m; i++) {
        s[i] = sin(MU_2PI * 0.05 * ts[i]);
    }

    matplotlibcpp::plot(ts, s);


    matplotlibcpp::show();




    return 0;
}