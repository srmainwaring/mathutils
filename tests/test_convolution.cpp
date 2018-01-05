//
// Created by frongere on 29/12/17.
//

#include "MathUtils.h"
#include "matplotlibcpp.h"

#include "boost/circular_buffer.hpp"

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
    auto dt = tr[1] - tr[0];
    auto r = Sinc(tr);

    // Creating the convolution engine
    Convolution<double> conv(dt, r);

    // Creating the signal into a circular buffer
    auto ts = linspace<double>(0., (double)m*100/(double)n, m);

    boost::circular_buffer<double> s(m, m, 0.); // Circular buffer initialized with capacity m and m values set to 0 (buffer is full)

//    std::vector<double> s(m);
    for (uint i=0; i<m; i++) {
        s.push_front(sin(MU_2PI * 0.05 * ts[i]));
    }

    // Vector construction
    std::vector<double> sv;
    for (unsigned int i=0; i<s.size(); i++) {
        sv.push_back(s[m-i-1]);
    }


//    matplotlibcpp::plot(ts, sv);
//    matplotlibcpp::show();




    return 0;
}