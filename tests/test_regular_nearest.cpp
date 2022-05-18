//
// Created by frongere on 18/05/22.
//

#include "MathUtils/RegularGridNearest.h"
#include "MathUtils/VectorGeneration.h"

#include <gtest/gtest.h>

using namespace mathutils;

TEST(RegularGridNearest, test_regular_grid_nearest) {

  using Array = boost::multi_array<int, 2>;
  Array var(boost::extents[3][3]);

  std::vector<double> xcoord = {0, 1, 2};
  std::vector<double> ycoord = {0, 1, 2};

  var[0][0] = 0;
  var[0][1] = 1;
  var[0][2] = 2;
  var[1][0] = 1;
  var[1][1] = 2;
  var[1][2] = 3;
  var[2][0] = 2;
  var[2][1] = 3;
  var[2][2] = 4;

  RegularGridNearest<int, 2> nearest;

  nearest.AddCoord(xcoord);
  nearest.AddCoord(ycoord);

  nearest.AddVar(var);

  ASSERT_EQ(nearest.Nearest({0., 0.}), 0);
  ASSERT_EQ(nearest.Nearest({0.25, 0.25}), 0);
  ASSERT_EQ(nearest.Nearest({0.5, 0.25}), 1);
  ASSERT_EQ(nearest.Nearest({0.75, 0.25}), 1);
  ASSERT_EQ(nearest.Nearest({1.25, 0.25}), 1);
  ASSERT_EQ(nearest.Nearest({1.75, 0.25}), 2);

  ASSERT_EQ(nearest.Nearest({0.25, 0.75}), 1);
  ASSERT_EQ(nearest.Nearest({0.75, 0.75}), 2);
  ASSERT_EQ(nearest.Nearest({1.25, 0.75}), 2);
  ASSERT_EQ(nearest.Nearest({1.75, 0.75}), 3);

  ASSERT_EQ(nearest.Nearest({0.25, 1.25}), 1);
  ASSERT_EQ(nearest.Nearest({0.75, 1.25}), 2);
  ASSERT_EQ(nearest.Nearest({1.25, 1.25}), 2);
  ASSERT_EQ(nearest.Nearest({1.75, 1.25}), 3);

  ASSERT_EQ(nearest.Nearest({0.25, 1.75}), 2);
  ASSERT_EQ(nearest.Nearest({0.75, 1.75}), 3);
  ASSERT_EQ(nearest.Nearest({1.25, 1.75}), 3);
  ASSERT_EQ(nearest.Nearest({1.75, 1.75}), 4);

//  ASSERT_EQ(nearest.Nearest({0.25, 0.25}), 0);
//  ASSERT_EQ(nearest.Nearest({0.25, 0.25}), 0);
//  ASSERT_EQ(nearest.Nearest({0.25, 0.25}), 0);
//  ASSERT_EQ(nearest.Nearest({0.25, 0.25}), 0);
//  ASSERT_EQ(nearest.Nearest({0.25, 0.25}), 0);
//  ASSERT_EQ(nearest.Nearest({0.25, 0.25}), 0);
//  ASSERT_EQ(nearest.Nearest({0.25, 0.25}), 0);
//  ASSERT_EQ(nearest.Nearest({0.25, 0.25}), 0);
//  ASSERT_EQ(nearest.Nearest({0.25, 0.25}), 0);
//  ASSERT_EQ(nearest.Nearest({0.25, 0.25}), 0);
//  ASSERT_EQ(nearest.Nearest({0.25, 0.25}), 0);
//  ASSERT_EQ(nearest.Nearest({0.25, 0.25}), 0);
//  ASSERT_EQ(nearest.Nearest({0.25, 0.25}), 0);
//  std::cout << nearest.Nearest({0.26, 0.54}) << std::endl; // 1
//  std::cout << nearest.Nearest({1.2, 1.9}) << std::endl;  // 3
//  std::cout << nearest.Nearest({1.49, 1.75}) << std::endl;  // 3
//  std::cout << nearest.Nearest({1.9, 1.75}) << std::endl;  // 4

}