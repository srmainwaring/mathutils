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

  // This test checks the computation of the double Chebyshev series.

  // Definition of a function to approximate by a double Chebyshev series.
  class FunctionTest : public Function2d<double> {
   public:
    /// This function evaluates the function at the point (x, y).
    double Evaluate(const double &x, const double &y) const override {
      return x * cos(y);
    }
    /// This function evaluates the x-derivative of the function at the point (x, y).
    double Evaluate_derivative_x(const double &x, const double &y) const {
      return cos(y);
    }
    /// This function evaluates the y-derivative of the function at the point (x, y).
    double Evaluate_derivative_y(const double &x, const double &y) const {
      return - x * sin(y);
    }
  };

  // Computation of the aij coefficients.
  FunctionTest myFunction;
  double xmin = -3.;
  double xmax = 2.;
  double ymin = -4.;
  double ymax = 8.;
  int order_x = 11;
  int order_y = 13;
  auto myDoubleChebyshevSeriesApprox = DoubleChebyshevSeriesApprox<double>(&myFunction, xmin, xmax, ymin, ymax, order_x, order_y);
  myDoubleChebyshevSeriesApprox.Computation_aij();

  // Tests.
  PrintHeader("Computation of the series approximation at a zero of the ChebyShev polynomial.");
  int r = 3;
  double tmp = cos(MU_PI_2 * (2. * r + 1.) / (order_x + 1.));
  double x = 0.5 * (xmax - xmin) * tmp + 0.5 * (xmax + xmin);
  int s = 3;
  tmp = cos(MU_PI_2 * (2. * s + 1.) / (order_y + 1.));
  double y = 0.5 * (ymax - ymin) * tmp + 0.5 * (ymax + ymin);
  std::cout << "Point: (x, y) = (" << x << ", " << y << ")" << std::endl;
  std::cout << "Analytical: " << myFunction.Evaluate(x, y) << std::endl;
  std::cout << "Double Chebyshev series approximation: " << myDoubleChebyshevSeriesApprox.Evaluate(x, y) << std::endl;
  assert(IsClose(myDoubleChebyshevSeriesApprox.Evaluate(x, y), myFunction.Evaluate(x, y)));

  PrintHeader("Computation of the series approximation at an arbitrary point.");
  std::cout << "Point: (1, 2)" << std::endl;
  std::cout << "Analytical: " << myFunction.Evaluate(1., 2.) << std::endl;
  std::cout << "Double Chebyshev series approximation: " << myDoubleChebyshevSeriesApprox.Evaluate(1., 2.) << std::endl;
  std::cout << "Relative error (%): " << 100 * abs((myDoubleChebyshevSeriesApprox.Evaluate(1., 2.) - myFunction.Evaluate(1., 2.))
  / myFunction.Evaluate(1., 2.)) << std::endl;

  PrintHeader("Computation of the x-derivative of the series approximation at an arbitrary point.");
  std::cout << "Point: (1, 2)" << std::endl;
  std::cout << "Analytical: " << myFunction.Evaluate_derivative_x(1., 2.) << std::endl;
  std::cout << "Double Chebyshev series approximation: " << myDoubleChebyshevSeriesApprox.Evaluate_derivative_x(1., 2.) << std::endl;
  std::cout << "Relative error (%): " << 100 * abs((myDoubleChebyshevSeriesApprox.Evaluate_derivative_x(1., 2.)
  - myFunction.Evaluate_derivative_x(1., 2.)) / myFunction.Evaluate_derivative_x(1., 2.)) << std::endl;

  PrintHeader("Computation of the y-derivative of the series approximation at an arbitrary point.");
  std::cout << "Point: (1, 2)" << std::endl;
  std::cout << "Analytical: " << myFunction.Evaluate_derivative_y(1., 2.) << std::endl;
  std::cout << "Double Chebyshev series approximation: " << myDoubleChebyshevSeriesApprox.Evaluate_derivative_y(1., 2.) << std::endl;
  std::cout << "Relative error (%): " << 100 * abs((myDoubleChebyshevSeriesApprox.Evaluate_derivative_y(1., 2.)
    - myFunction.Evaluate_derivative_y(1., 2.)) / myFunction.Evaluate_derivative_y(1., 2.)) << std::endl;

}