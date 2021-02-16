//
// Created by pierre-yves on 01/12/2020.
//

#include <iostream>
#include "MathUtils/MathUtils.h"

using namespace mathutils;

void PrintHeader(std::string title) {
  std::cout << "\n=====================================================================" << std::endl;
  std::cout << "    " << title << std::endl;
  std::cout << "=====================================================================" << std::endl;
}

void PrintInfo(std::string info) {
  std::cout << info << ":" << std::endl;
}

int main(int argc, char* argv[]) {

  // This test checks the computation of the double and triple Chebyshev series.

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //                                              DOUBLE - Closed segments.
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // Definition of a function to approximate by a double Chebyshev series.
  class Function2dTest : public Function2d<double> {
   public:
    /// This method evaluates the function at the point (x, y).
    double Evaluate(const double &x, const double &y) const override {
      return x * cos(y);
    }
    /// This method evaluates the x-derivative of the function at the point (x, y).
    double Evaluate_derivative_x(const double &x, const double &y) const {
      return cos(y);
    }
    /// This method evaluates the y-derivative of the function at the point (x, y).
    double Evaluate_derivative_y(const double &x, const double &y) const {
      return - x * sin(y);
    }
  };

  // Computation of the aij coefficients.
  Function2dTest myFunction2d;
  double xmin = -3.;
  double xmax = 2.;
  double ymin = -4.;
  double ymax = 8.;
  int order_x = 11;
  int order_y = 13;
  auto myDoubleChebyshevSeriesApprox = DoubleChebyshevSeriesApproxClosedSegment<double>(&myFunction2d, xmin, xmax, ymin, ymax, order_x, order_y);
  myDoubleChebyshevSeriesApprox.Computation_aij();

  // Tests.
  PrintHeader("Computation of the DOUBLE series approximation at a zero of the ChebyShev polynomial.");
  int r = 3;
  double tmp = cos(MU_PI_2 * (2. * r + 1.) / (order_x + 1.));
  double x = 0.5 * (xmax - xmin) * tmp + 0.5 * (xmax + xmin);
  int s = 3;
  tmp = cos(MU_PI_2 * (2. * s + 1.) / (order_y + 1.));
  double y = 0.5 * (ymax - ymin) * tmp + 0.5 * (ymax + ymin);
  std::cout << "Point: (x, y) = (" << x << ", " << y << ")" << std::endl;
  std::cout << "Analytical: " << myFunction2d.Evaluate(x, y) << std::endl;
  std::cout << "Double Chebyshev series approximation: " << myDoubleChebyshevSeriesApprox.Evaluate(x, y) << std::endl;
  assert(IsClose(myDoubleChebyshevSeriesApprox.Evaluate(x, y), myFunction2d.Evaluate(x, y)));

  PrintHeader("Computation of the DOUBLE series approximation at an arbitrary point.");
  std::cout << "Point: (1, 2)" << std::endl;
  std::cout << "Analytical: " << myFunction2d.Evaluate(1., 2.) << std::endl;
  std::cout << "Double Chebyshev series approximation: " << myDoubleChebyshevSeriesApprox.Evaluate(1., 2.) << std::endl;
  std::cout << "Relative error (%): " << 100 * abs((myDoubleChebyshevSeriesApprox.Evaluate(1., 2.) - myFunction2d.Evaluate(1., 2.))
  / myFunction2d.Evaluate(1., 2.)) << std::endl;
  assert(IsClose(myDoubleChebyshevSeriesApprox.Evaluate(1., 2.), myFunction2d.Evaluate(1., 2.), 10e-3));

  PrintHeader("Computation of the x-derivative of the DOUBLE series approximation at an arbitrary point.");
  std::cout << "Point: (1, 2)" << std::endl;
  std::cout << "Analytical: " << myFunction2d.Evaluate_derivative_x(1., 2.) << std::endl;
  std::cout << "Double Chebyshev series approximation: " << myDoubleChebyshevSeriesApprox.Evaluate_derivative_x(1., 2.) << std::endl;
  std::cout << "Relative error (%): " << 100 * abs((myDoubleChebyshevSeriesApprox.Evaluate_derivative_x(1., 2.)
  - myFunction2d.Evaluate_derivative_x(1., 2.)) / myFunction2d.Evaluate_derivative_x(1., 2.)) << std::endl;
  assert(IsClose(myDoubleChebyshevSeriesApprox.Evaluate_derivative_x(1., 2.), myFunction2d.Evaluate_derivative_x(1., 2.), 10e-3));

  PrintHeader("Computation of the y-derivative of the DOUBLE series approximation at an arbitrary point.");
  std::cout << "Point: (1, 2)" << std::endl;
  std::cout << "Analytical: " << myFunction2d.Evaluate_derivative_y(1., 2.) << std::endl;
  std::cout << "Double Chebyshev series approximation: " << myDoubleChebyshevSeriesApprox.Evaluate_derivative_y(1., 2.) << std::endl;
  std::cout << "Relative error (%): " << 100 * abs((myDoubleChebyshevSeriesApprox.Evaluate_derivative_y(1., 2.)
    - myFunction2d.Evaluate_derivative_y(1., 2.)) / myFunction2d.Evaluate_derivative_y(1., 2.)) << std::endl;
  assert(IsClose(myDoubleChebyshevSeriesApprox.Evaluate_derivative_y(1., 2.), myFunction2d.Evaluate_derivative_y(1., 2.), 10e-3));

  std::cout << "\n\n" << std::endl;

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //                                          DOUBLE - Half-open segments
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // Definition of a second function to approximate by a double Chebyshev series.
  class Function2dTestHalfOpen : public Function2d<double> {
   public:
    /// This method evaluates the function at the point (x, y).
    double Evaluate(const double &x, const double &y) const override {
      return 1. / (2. * x + 3. * y);
    }
    /// This method evaluates the x-derivative of the function at the point (x, y).
    double Evaluate_derivative_x(const double &x, const double &y) const {
      return - 2. / pow(2. * x + 3. * y, 2.);
    }
    /// This method evaluates the y-derivative of the function at the point (x, y).
    double Evaluate_derivative_y(const double &x, const double &y) const {
      return - 3. / pow(2. * x + 3. * y, 2.);
    }
  };

  // Computation of the aij coefficients.
  Function2dTestHalfOpen myFunction2dHalfOpen;
  xmin = 1.;
  ymin = 1.;
  order_x = 50;
  order_y = 50;
  x = 150;
  y = 150;

  // Tests.
  auto myDoubleChebyshevSeriesApproxHalfOpen = DoubleChebyshevSeriesApproxHalfOpenSegment<double>(&myFunction2dHalfOpen, xmin, ymin, order_x, order_y);
  myDoubleChebyshevSeriesApproxHalfOpen.Computation_aij();
  PrintHeader("Computation of the DOUBLE series approximation at an arbitrary point.");
  std::cout << "Point: (150, 150)" << std::endl;

  std::cout << "Analytical: " << myFunction2dHalfOpen.Evaluate(x, y) << std::endl;
  std::cout << "Double Chebyshev series approximation: " << myDoubleChebyshevSeriesApproxHalfOpen.Evaluate(x, y) << std::endl;
  std::cout << "Relative error (%): " << 100 * abs((myDoubleChebyshevSeriesApproxHalfOpen.Evaluate(x, y) - myFunction2dHalfOpen.Evaluate(x, y))
                                                   / myFunction2dHalfOpen.Evaluate(x, y)) << std::endl;
  assert(IsClose(myDoubleChebyshevSeriesApproxHalfOpen.Evaluate(x, y), myFunction2dHalfOpen.Evaluate(x, y), 10e-3));

  PrintHeader("Computation of the x-derivative of the DOUBLE series approximation at an arbitrary point.");
  std::cout << "Point: (150, 150)" << std::endl;

  std::cout << "Analytical: " << myFunction2dHalfOpen.Evaluate_derivative_x(x, y) << std::endl;
  std::cout << "Double Chebyshev series approximation: " << myDoubleChebyshevSeriesApproxHalfOpen.Evaluate_derivative_x(x, y) << std::endl;
  std::cout << "Relative error (%): " << 100 * abs((myDoubleChebyshevSeriesApproxHalfOpen.Evaluate_derivative_x(x, y) - myFunction2dHalfOpen.Evaluate_derivative_x(x, y))
                                                   / myFunction2dHalfOpen.Evaluate_derivative_x(x, y)) << std::endl;
  assert(IsClose(myDoubleChebyshevSeriesApproxHalfOpen.Evaluate_derivative_x(x, y), myFunction2dHalfOpen.Evaluate_derivative_x(x, y), 10e-3));

  PrintHeader("Computation of the y-derivative of the DOUBLE series approximation at an arbitrary point.");
  std::cout << "Point: (150, 150)" << std::endl;

  std::cout << "Analytical: " << myFunction2dHalfOpen.Evaluate_derivative_y(x, y) << std::endl;
  std::cout << "Double Chebyshev series approximation: " << myDoubleChebyshevSeriesApproxHalfOpen.Evaluate_derivative_y(x, y) << std::endl;
  std::cout << "Relative error (%): " << 100 * abs((myDoubleChebyshevSeriesApproxHalfOpen.Evaluate_derivative_y(x, y) - myFunction2dHalfOpen.Evaluate_derivative_y(x, y))
                                                   / myFunction2dHalfOpen.Evaluate_derivative_y(x, y)) << std::endl;
  assert(IsClose(myDoubleChebyshevSeriesApproxHalfOpen.Evaluate_derivative_y(x, y), myFunction2dHalfOpen.Evaluate_derivative_y(x, y), 10e-3));

  std::cout << "\n\n" << std::endl;

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //                                              TRIPLE
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // Definition of a function to approximate by a double Chebyshev series.
  class Function3dTest : public Function3d<double> {
   public:
    /// This method evaluates the function at the point (x, y, z).
    double Evaluate(const double &x, const double &y, const double &z) const override {
      return x * cos(y) * sin(z);
    }
    /// This method evaluates the x-derivative of the function at the point (x, y, z).
    double Evaluate_derivative_x(const double &x, const double &y, const double &z) const {
      return cos(y) * sin(z);
    }
    /// This method evaluates the y-derivative of the function at the point (x, y, z).
    double Evaluate_derivative_y(const double &x, const double &y, const double &z) const {
      return - x * sin(y) * sin(z);
    }
    /// This method evaluates the z-derivative of the function at the point (x, y, z).
    double Evaluate_derivative_z(const double &x, const double &y, const double &z) const {
      return x * cos(y) * cos(z);
    }
  };

  // Computation of the aij coefficients.
  Function3dTest myFunction3d;
  xmin = -3.;
  xmax = 2.;
  ymin = -4.;
  ymax = 8.;
  double zmin = -5.;
  double zmax = 4.;
  order_x = 11;
  order_y = 13;
  int order_z = 10;
  auto myTripleChebyshevSeriesApprox = TripleChebyshevSeriesApprox<double>(&myFunction3d, xmin, xmax, ymin, ymax, zmin,
                                                                           zmax, order_x, order_y, order_z);
  myTripleChebyshevSeriesApprox.Computation_aijk();

  // Tests.
  PrintHeader("Computation of the TRIPLE series approximation at a zero of the ChebyShev polynomial.");
  r = 3;
  tmp = cos(MU_PI_2 * (2. * r + 1.) / (order_x + 1.));
  x = 0.5 * (xmax - xmin) * tmp + 0.5 * (xmax + xmin);
  s = 3;
  tmp = cos(MU_PI_2 * (2. * s + 1.) / (order_y + 1.));
  y = 0.5 * (ymax - ymin) * tmp + 0.5 * (ymax + ymin);
  int t = 3;
  tmp = cos(MU_PI_2 * (2. * t + 1.) / (order_z + 1.));
  double z = 0.5 * (zmax - zmin) * tmp + 0.5 * (zmax + zmin);
  std::cout << "Point: (x, y, z) = (" << x << ", " << y << ", " << z << ")" << std::endl;
  std::cout << "Analytical: " << myFunction3d.Evaluate(x, y, z) << std::endl;
  std::cout << "Triple Chebyshev series approximation: " << myTripleChebyshevSeriesApprox.Evaluate(x, y, z) << std::endl;
  assert(IsClose(myTripleChebyshevSeriesApprox.Evaluate(x, y, z), myFunction3d.Evaluate(x, y, z)));

  PrintHeader("Computation of the TRIPLE series approximation at an arbitrary point.");
  std::cout << "Point: (1, 2, 3)" << std::endl;
  std::cout << "Analytical: " << myFunction3d.Evaluate(1., 2., 3) << std::endl;
  std::cout << "Tripe Chebyshev series approximation: " << myTripleChebyshevSeriesApprox.Evaluate(1., 2., 3.) << std::endl;
  std::cout << "Relative error (%): " << 100 * abs((myTripleChebyshevSeriesApprox.Evaluate(1., 2., 3.) - myFunction3d.Evaluate(1., 2., 3.))
                                                   / myFunction3d.Evaluate(1., 2., 3.)) << std::endl;

  PrintHeader("Computation of the x-derivative of the TRIPLE series approximation at an arbitrary point.");
  std::cout << "Point: (1, 2, 3)" << std::endl;
  std::cout << "Analytical: " << myFunction3d.Evaluate_derivative_x(1., 2., 3.) << std::endl;
  std::cout << "Triple Chebyshev series approximation: " << myTripleChebyshevSeriesApprox.Evaluate_derivative_x(1., 2., 3.) << std::endl;
  std::cout << "Relative error (%): " << 100 * abs((myTripleChebyshevSeriesApprox.Evaluate_derivative_x(1., 2., 3.)
                                                    - myFunction3d.Evaluate_derivative_x(1., 2., 3.)) / myFunction3d.Evaluate_derivative_x(1., 2., 3.)) << std::endl;

  PrintHeader("Computation of the y-derivative of the TRIPLE series approximation at an arbitrary point.");
  std::cout << "Point: (1, 2)" << std::endl;
  std::cout << "Analytical: " << myFunction3d.Evaluate_derivative_y(1., 2., 3.) << std::endl;
  std::cout << "Triple Chebyshev series approximation: " << myTripleChebyshevSeriesApprox.Evaluate_derivative_y(1., 2., 3.) << std::endl;
  std::cout << "Relative error (%): " << 100 * abs((myTripleChebyshevSeriesApprox.Evaluate_derivative_y(1., 2., 3.)
                                                    - myFunction3d.Evaluate_derivative_y(1., 2., 3.)) / myFunction3d.Evaluate_derivative_y(1., 2., 3.)) << std::endl;

  PrintHeader("Computation of the z-derivative of the TRIPLE series approximation at an arbitrary point.");
  std::cout << "Point: (1, 2, 3)" << std::endl;
  std::cout << "Analytical: " << myFunction3d.Evaluate_derivative_z(1., 2., 3.) << std::endl;
  std::cout << "Triple Chebyshev series approximation: " << myTripleChebyshevSeriesApprox.Evaluate_derivative_z(1., 2., 3.) << std::endl;
  std::cout << "Relative error (%): " << 100 * abs((myTripleChebyshevSeriesApprox.Evaluate_derivative_z(1., 2., 3.)
                                                    - myFunction3d.Evaluate_derivative_z(1., 2., 3.)) / myFunction3d.Evaluate_derivative_z(1., 2., 3.)) << std::endl;


}