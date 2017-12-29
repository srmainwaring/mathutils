//
// Created by frongere on 29/12/17.
//

#include "MathUtils.h"
#include "matplotlibcpp.h"

using namespace mathutils;


int main(int argc, char const* argv[])
{

    int n = 100;
    auto x = linspace<double>(0, 100, n);

//    auto y = Door<double>(x, 25, 40, 5.5);
//    auto y = Dirac<double>(x, 0);

    auto y = Sinc(x);

    Spline<double, 1, 3> s(x, y);


    int ni = 1000;
    auto xi = linspace<double>(0, 100, ni);
    std::vector<double> yi(ni);
    for (uint i=0; i<ni; i++) {
        yi[i] = s(xi[i]);
    }


    matplotlibcpp::plot(x, y);
    matplotlibcpp::plot(xi, yi);
    matplotlibcpp::show();
}