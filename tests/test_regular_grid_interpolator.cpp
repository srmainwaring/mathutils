//
// Created by frongere on 27/12/2021.
//

#include "MathUtils/RegularGridInterpolator.h"
#include "MathUtils/VectorGeneration.h"

#include <gtest/gtest.h>

using namespace mathutils;

inline double func(const double &x,
                   const double &y,
                   const double &z) {
  return x + 2 * std::pow(y, 2) + 3 * std::pow(z, 3);
}

double frand(double f_min, double f_max) {
  double f = (double) rand() / RAND_MAX;
  return f_min + f * (f_max - f_min);
}


TEST(RegularGridInterpolator, test_regular_grid_interpolator) {

  using Array3D = boost::multi_array<double, 3>;
  using index = Array3D::index;

  int nx = 100;
  int ny = 200;
  int nz = 300;

  double xmin = -1;
  double xmax = 1;
  double ymin = -2;
  double ymax = 2;
  double zmin = -3;
  double zmax = 3;

  // Building dataset to interpolate
  auto xcoord = linspace<double>(xmin, xmax, nx);
  auto ycoord = linspace<double>(ymin, ymax, ny);
  auto zcoord = linspace<double>(zmin, zmax, nz);

  Array3D var(boost::extents[nx][ny][nz]);

  for (index ix = 0; ix < nx; ++ix) {
    for (index iy = 0; iy < ny; ++iy) {
      for (index iz = 0; iz < nz; ++iz) {
        var[ix][iy][iz] = func(xcoord[ix], ycoord[iy], zcoord[iz]);
      }
    }
  }

  // Building interpolator
  RegularGridInterpolator<double, 3> interpolator;

  interpolator.AddCoord(xcoord);
  interpolator.AddCoord(ycoord);
  interpolator.AddCoord(zcoord);

  interpolator.AddVar(var);

  int ntest = 1000000;
  for (int i = 0; i < ntest; ++i) {
    double x = frand(xmin, xmax);
    double y = frand(ymin, ymax);
    double z = frand(zmin, zmax);

    ASSERT_NEAR(interpolator.Interp({x, y, z}), func(x, y, z), 1e-1);

  }
}

//TEST(RegularGridInterpolator, test_regular_grid_nearest_double) {
//
//  using Array = boost::multi_array<double, 2>;
//  Array var(boost::extents[3][3]);
//
//  std::vector<double> xcoord = {0, 1, 2};
//  std::vector<double> ycoord = {0, 1, 2};
//
//  var[0][0] = 0;
//  var[0][1] = 1;
//  var[0][2] = 2;
//  var[1][0] = 1;
//  var[1][1] = 2;
//  var[1][2] = 3;
//  var[2][0] = 2;
//  var[2][1] = 3;
//  var[2][2] = 4;
//
//  RegularGridInterpolator<double, 2> interpolator;
//
//  interpolator.AddCoord(xcoord);
//  interpolator.AddCoord(ycoord);
//
//  interpolator.AddVar(var);
//
//  std::cout << interpolator.Interp({0.26, 0.25}) << std::endl;
//  std::cout << interpolator.Interp({0.9, 1.9}) << std::endl;
////  std::cout << interpolator.Nearest({0.26, 0.25}) << std::endl;
////  std::cout << interpolator.Nearest({0.9, 1.9}) << std::endl;
//
//
//  // Attention, si on depasse les ranges de l'interpolateur, on est clamp et ça dit rien de special...
//
////  std::cout << interpolator.Nearest({0.5, 0.5}) << std::endl;
//
//
//}
