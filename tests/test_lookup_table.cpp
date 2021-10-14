//
// Created by frongere on 16/11/17.
//


//#include <iostream>
#include "MathUtils/MathUtils.h"

#include <gtest/gtest.h>


using namespace mathutils;


TEST(LookupTable, OldMain) {
    auto x = linspace<double>(0, 100, 100);
    auto cx = linspace<double>(10, 90, 100);
    auto cy = linspace<double>(0.2, 87., 100);

    auto lookup = LookupTable1D<double>(LINEAR);
    lookup.SetX(x);
    lookup.AddY("cx", cx);
    lookup.AddY("cy", cy);

    auto xInterp = linspace<double>(0.1, 99.9, 200);

    lookup.Eval("cx", xInterp);

    for (auto xi : xInterp) {

        ASSERT_NEAR(10 + 80 * xi/100, lookup.Eval("cx", xi), 1e-8);
        ASSERT_NEAR(0.2 + 86.8 * xi/100, lookup.Eval("cy", xi), 1e-8);

    }
}
