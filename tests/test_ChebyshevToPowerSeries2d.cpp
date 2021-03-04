//
// Created by pierre-yves on 04/03/2021.
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

  // This test checks the transformation from double Chebyshev series to double power series.

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //                                              Closed segments.
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // Definition of a function to approximate by a double Chebyshev series.
  class Function2dTest : public Function2d<double> {
   public:
    /// This method evaluates the function at the point (x, y).
    double Evaluate(const double &x, const double &y) const override {
      return 1. / (2. * x + 3. * y);
    }

    /// This method evaluates the x-derivative of the function at the point (x, y).
    double Evaluate_derivative_x(const double &x, const double &y) const {
      return -2. / pow(2. * x + 3. * y, 2.);
    }

    /// This method evaluates the y-derivative of the function at the point (x, y).
    double Evaluate_derivative_y(const double &x, const double &y) const {
      return -3. / pow(2. * x + 3. * y, 2.);
    }
  };

  // aij coefficients.
  Function2dTest myFunction2d;
  double xmin = 1.;
  double xmax = 10.;
  double ymin = 1.;
  double ymax = 10.;
  int order_x = 11;
  int order_y = 13;
  auto myChebyshevSeries2dClosed = std::make_shared<ChebyshevSeries2dClosed<double>>(&myFunction2d, xmin, xmax, ymin, ymax, order_x, order_y);
  myChebyshevSeries2dClosed->Compute_aij();

  // From Chebyshev series to power series.
  auto myPowerSeries2dClosed = ChebyshevToPowerSeries2dClosed<double>(myChebyshevSeries2dClosed);

  PrintHeader("Closed segments - Function.");
  double x = 3.;
  double y = 2.;
  std::cout << "Point: (x, y) = (" << x << ", " << y << ")" << std::endl;
  double ana = myFunction2d.Evaluate(x, y);
  std::cout << "Analytical: " << ana << std::endl;
  double chebyshev = myChebyshevSeries2dClosed->Evaluate(x, y);
  std::cout << "Double Chebyshev series approximation: " << chebyshev << std::endl;
  double power = myPowerSeries2dClosed.Evaluate(x, y);
  std::cout << "Double Power series approximation: " << power << std::endl;
  std::cout << "Relative error (%): " << 100 * abs((chebyshev - power) / power) << std::endl;
  assert(IsClose(ana, power, 10e-3));
  assert(IsClose(chebyshev, power, 10e-3));

  PrintHeader("Closed segments - x derivative.");
  x = 3.;
  y = 2.;
  std::cout << "Point: (x, y) = (" << x << ", " << y << ")" << std::endl;
  ana = myFunction2d.Evaluate_derivative_x(x, y);
  std::cout << "Analytical: " << ana << std::endl;
  chebyshev = myChebyshevSeries2dClosed->Evaluate_derivative_x(x, y);
  std::cout << "Double Chebyshev series approximation: " << chebyshev << std::endl;
  power = myPowerSeries2dClosed.Evaluate_derivate_x(x, y);
  std::cout << "Double Power series approximation: " << power << std::endl;
  std::cout << "Relative error (%): " << 100 * abs((chebyshev - power) / power) << std::endl;
  assert(IsClose(ana, power, 10e-3));
  assert(IsClose(chebyshev, power, 10e-3));

  PrintHeader("Closed segments - y derivative.");
  x = 3.;
  y = 2.;
  std::cout << "Point: (x, y) = (" << x << ", " << y << ")" << std::endl;
  ana = myFunction2d.Evaluate_derivative_y(x, y);
  std::cout << "Analytical: " << ana << std::endl;
  chebyshev = myChebyshevSeries2dClosed->Evaluate_derivative_y(x, y);
  std::cout << "Double Chebyshev series approximation: " << chebyshev << std::endl;
  power = myPowerSeries2dClosed.Evaluate_derivate_y(x, y);
  std::cout << "Double Power series approximation: " << power << std::endl;
  std::cout << "Relative error (%): " << 100 * abs((chebyshev - power) / power) << std::endl;
  assert(IsClose(ana, power, 10e-3));
  assert(IsClose(chebyshev, power, 10e-3));

}
