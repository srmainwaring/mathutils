//
// Created by frongere on 18/12/17.
//

#include "MathUtils/MathUtils.h"

using namespace mathutils;


// Defining a function to integrate
double myFunction(double x) {
    if (x == 0.) {
        return 1.;
    } else {
        return sin(x) / x;
    }
}


int main(int argc, char* argv[]) {

    // Defining the integration.
    // Here, by default we use the trapezoidal method
    auto myIntegrator = Integrate1d<double>(myFunction, 0, 10, 100);
    std::cout << myIntegrator.Get() << std::endl;


    // Defining a vector
    auto x = linspace<double>(0, 10, 100);
    std::vector<double> y;
    for (auto val : x) {
        y.push_back(myFunction(val));
    }

    auto myIntegratorVect = Integrate1d<double>(y, 0, 10, 100);
    std::cout << myIntegratorVect.Get() << std::endl;








    return 0;
}