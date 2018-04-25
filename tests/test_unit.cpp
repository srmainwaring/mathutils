//
// Created by frongere on 16/11/17.
//

#include <iostream>
#include "MathUtils/MathUtils.h"

using namespace mathutils;

int main(int argc, char* argv[]) {

    assert(deg2rad(180.) == M_PI);
    assert(rad2deg<double>(M_PI) == 180);

    // Pas tres utilse de tout faire...

    return 0;
}