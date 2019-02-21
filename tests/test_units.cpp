//
// Created by frongere on 16/11/18.
//

#include "units.h"

#include "MathUtils/MathUtils.h"

using namespace mathutils;

using namespace units::literals;

using namespace units;
using namespace units::length;
using namespace units::time;
using namespace units::area;
using namespace units::velocity;
using namespace units::angle;
using namespace units::angular_velocity;


int main() {

    // Distances
    meter_t d = 2._m;

    std::cout << (double)d << std::endl;


    std::cout << d * 2. << std::endl;

    double dd = 5.;
    std::cout << (meter_t)dd << std::endl;



//    Vector2d<meter_t> pos;
//    pos << 1_m, 2_m;
//    pos[0] = 1_m;




//    Vector3d<double> pos;
//    pos << 1, 2;

//    std::cout << pos << std::endl;



    // Essai avec des vitesses
    meters_per_second_t v = 2.0_mps;

    kilometers_per_hour_t vkmh = v;
    std::cout << vkmh << std::endl;


    // angles



}