//
// Created by frongere on 16/11/17.
//

#include <iostream>
#include "MathUtils/MathUtils.h"
#include "gtest/gtest.h"

#define N 100

using namespace mathutils;

TEST(Interp1d, Double) {

  // Building the x coords as a shared pointer
  auto x = std::make_shared<std::vector<double>>(
      linspace(M_PI_2, 4 * M_PI, N - 1)
  );

  // Building the data
  auto y = std::make_shared<std::vector<double>>();
  y->reserve(x->size());
  double val;
  for (unsigned long i = 0; i < x->size(); i++) {
    val = sin(x->at(i));
    y->push_back(val);
  }

  // Create the linear interpolator
  Interp1dLinear<double, double> interpolator;
  interpolator.Initialize(x, y);

  // Test of the Eval method for one scalar
  EXPECT_EQ(interpolator.Eval(5.3333), interpolator(5.3333));
  EXPECT_NEAR(sin(5.3333), interpolator(5.3333), 1E-2);

  // Test exit if outside of ranges
  EXPECT_EXIT(interpolator(0.), ::testing::ExitedWithCode(1), ".*");

  // Test for a vector of x coords
  auto x_interp = linspace(M_PI_2, 4 * M_PI, 1000 * N);
  // Using only the overloaded call operator for vector values
  auto y_interp = interpolator(x_interp);

  // Interpolator saturates outside of the ranges
  Interp1dLinearSaturate<double, double> saturator;
  saturator.Initialize(x, y);

  EXPECT_EQ(saturator(MU_PI_2), saturator(0.));
  EXPECT_EQ(saturator(4 * MU_PI), saturator(5 * MU_PI));

  // Linear extrapolator
  Interp1dLinearExtrapolate<double, double> extrapolator;
  extrapolator.Initialize(x, y);

  auto a = (sin(x->at(1)) - sin(x->at(0))) / (x->at(1) - x->at(0));
  auto b = sin(x->at(0)) - a * x->at(0);

  EXPECT_NEAR(b, extrapolator(0.), 1E-8);

}

TEST(Interp1d, Complex) {

  // Building the x coords as a shared pointer
  auto x = std::make_shared<std::vector<double>>(
      linspace(M_PI, 4 * M_PI, N - 1)
  );

  // Building the data
  auto y = std::make_shared<std::vector<std::complex<double>>>();
  y->reserve(x->size());
  std::complex<double> val;
  for (unsigned long i = 0; i < x->size(); i++) {
    val = exp(MU_JJ * x->at(i));
    y->push_back(val);
  }

  // Create the interpolation
  Interp1dLinear<double, std::complex<double>> interpolator;

  interpolator.Initialize(x, y);

  // Test of the Eval method for one scalar
  EXPECT_EQ(interpolator.Eval(5.3333), interpolator(5.3333));
  EXPECT_NEAR(exp(MU_JJ * 5.3333).real(), interpolator(5.3333).real(), 1E-2);
  EXPECT_NEAR(exp(MU_JJ * 5.3333).imag(), interpolator(5.3333).imag(), 1E-2);

  // Test for a vector of x coords
  auto x_interp = linspace(M_PI, 4 * M_PI, 1000 * N);
  // Using only the overloaded call operator for vector values
  auto y_interp = interpolator(x_interp);
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}