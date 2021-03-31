//
// Created by frongere on 29/12/17.
//

#include "MathUtils/MathUtils.h"
//#include "matplotlibcpp.h"

using namespace mathutils;


int main(int argc, char const* argv[])
{

    int n = 1000;
    auto x = linspace<double>(0, 100, n);

    auto y = Heaviside<double>(x, 20);
//    auto y = Door<double>(x, 25, 40, 5.5);
//    auto y = Dirac<double>(x, 10);
//
//    auto y = Sinc(x);

    Spline<double, 1, 3> s(x, y);


    // TODO: ajouter methode pour l'evaluation de la spline sur un echantillopn std::vector...
    int ni = 10000;
    auto xi = linspace<double>(0, 100, ni);
    std::vector<double> yi(ni);
    for (unsigned int i=0; i<ni; i++) {
        yi[i] = s(xi[i]);
    }


//    matplotlibcpp::plot(x, y);
//    matplotlibcpp::plot(xi, yi);
//    matplotlibcpp::show();
}
