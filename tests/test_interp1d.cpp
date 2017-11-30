//
// Created by frongere on 16/11/17.
//

#include <iostream>
#include "MathUtils.h"

#define N 100

using namespace mathutils;


void testInterpDouble() {
    // Building the x coords as a shared pointer
    auto x = std::make_shared<std::vector<double>>(
            linspace(M_PI, 4*M_PI, N-1)
    );

    // Building the data
    auto y = std::make_shared<std::vector<double>>();
    y->reserve(x->size());
    double val;
    for (unsigned long i = 0; i < x->size(); i++) {
        val = sin(x->at(i));
        y->push_back( val );
    }

    // Create the interpolation
    Interp1dLinear<double, double> interpolator;

    interpolator.Initialize(x, y);

    // Test of the Eval method for one scalar
    auto y0 = interpolator.Eval(5.3333);
    // Test of the call operator for one scalar
    auto y1 = interpolator(5.3333);

    assert(IsClose(y0, y1));
//    assert(IsClose(y0, -0.8133409832926298));

    // Test for a vector of x coords
    auto x_interp = linspace(M_PI, 4*M_PI, 1000*N);
    // Using only the overloaded call operator for vector values
    auto y_interp = interpolator(x_interp);
}

void testInterpComplex() {
    // Building the x coords as a shared pointer
    auto x = std::make_shared<std::vector<double>>(
            linspace(M_PI, 4*M_PI, N-1)
    );

    // Building the data
    auto y = std::make_shared<std::vector<std::complex<double>>>();
    y->reserve(x->size());
    std::complex<double> val;
    for (unsigned long i = 0; i < x->size(); i++) {
        val = exp(MU_JJ * x->at(i));
        y->push_back( val );
    }

    // Create the interpolation
    Interp1dLinear<double, std::complex<double>> interpolator;

    interpolator.Initialize(x, y);

    // Test of the Eval method for one scalar
    auto y0 = interpolator.Eval(5.3333);
    // Test of the call operator for one scalar
    auto y1 = interpolator(5.3333);

//    assert(IsClose(y0, y1));
//    assert(IsClose(y0, -0.8133409832926298));

    // Test for a vector of x coords
    auto x_interp = linspace(M_PI, 4*M_PI, 1000*N);
    // Using only the overloaded call operator for vector values
    auto y_interp = interpolator(x_interp);
}


int main(int argc, char* argv[]) {

    testInterpDouble();
    testInterpComplex();

    return 0;
}