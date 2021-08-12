//
// Created by frongere on 16/11/17.
//


//#include <iostream>
#include "MathUtils/MathUtils.h"

using namespace mathutils;

int main(int argc, char* argv[]) {

    auto x = linspace<double>(0, 100, 100);
    auto cx = linspace<double>(10, 90, 100);
    auto cy = linspace<double>(0.2, 87., 100);

    auto lookup = LookupTable1D<double>();
    lookup.SetX(x);
    lookup.AddY("cx", cx);
    lookup.AddY("cy", cy);

    auto xInterp = linspace<double>(0.1, 99.9, 200);

    lookup.Eval("cx", xInterp);

    for (auto xi : xInterp) {

        lookup.Eval("cx", xi);
        lookup.Eval("cy", xi);

    }


    return 0;
}
