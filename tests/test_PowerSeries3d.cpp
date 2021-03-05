//
// Created by pierre-yves on 05/03/2021.
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

  // This test checks the transformation from triple Chebyshev series to triple power series.

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //                                              Closed segments
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // Definition of a function to approximate by a double Chebyshev series.
  class Function3dTest : public Function3d<double> {
   public:
    /// This method evaluates the function at the point (x, y, z).
    double Evaluate(const double &x, const double &y, const double &z) const override {
      return 1. / (2. * x + 3. * y + 4. * z);
    }

    /// This method evaluates the x-derivative of the function at the point (x, y, z).
    double Evaluate_derivative_x(const double &x, const double &y, const double &z) const {
      return -2. / pow(2. * x + 3. * y + 4. * z, 2.);
    }

    /// This method evaluates the y-derivative of the function at the point (x, y, z).
    double Evaluate_derivative_y(const double &x, const double &y, const double &z) const {
      return -3. / pow(2. * x + 3. * y + 4. * z, 2.);
    }

    /// This method evaluates the z-derivative of the function at the point (x, y, z).
    double Evaluate_derivative_z(const double &x, const double &y, const double &z) const {
      return -4. / pow(2. * x + 3. * y + 4. * z, 2.);
    }
  };

  // aij coefficients.
  Function3dTest myFunction3d;
  double xmin = 1.;
  double xmax = 6.;
  double ymin = 1.;
  double ymax = 6.;
  double zmin = 1.;
  double zmax = 6.;
  int order_x = 11;
  int order_y = 13;
  int order_z = 10;
  auto myChebyshevSeries3dClosed = ChebyshevSeries3dClosed<double>(&myFunction3d, xmin, xmax, ymin, ymax, zmin, zmax,
                                                                   order_x, order_y, order_z);
  myChebyshevSeries3dClosed.Compute_aijk();

  // From Chebyshev series to power series.
  auto bijk = myChebyshevSeries3dClosed.Compute_bijk();
  auto myPowerSeries3dClosed = PowerSeries3dClosed<double>(bijk, xmin, xmax, ymin, ymax, zmin, zmax);

  PrintHeader("Closed segments - Function.");
  double x = 3.;
  double y = 2.;
  double z = 5.;
  std::cout << "Point: (x, y, z) = (" << x << ", " << y << ", " << z << ")" << std::endl;
  double ana = myFunction3d.Evaluate(x, y, z);
  std::cout << "Analytical: " << ana << std::endl;
  double chebyshev = myChebyshevSeries3dClosed.Evaluate(x, y, z);
  std::cout << "Tripe Chebyshev series approximation: " << chebyshev << std::endl;
  std::cout << "Relative error (%): " << 100 * abs((chebyshev - ana) / ana) << std::endl;
  double power = myPowerSeries3dClosed.Evaluate(x, y, z);
  std::cout << "Triple Power series approximation: " << power << std::endl;
  std::cout << "Relative error (%): " << 100 * abs((power - ana) / ana) << std::endl;
  assert(IsClose(ana, chebyshev, 10e-3));
  assert(IsClose(ana, power, 10e-3));

  PrintHeader("Closed segments - x derivative.");
  ana = myFunction3d.Evaluate_derivative_x(x, y, z);
  std::cout << "Analytical: " << ana << std::endl;
  chebyshev = myChebyshevSeries3dClosed.Evaluate_derivative_x(x, y, z);
  std::cout << "Tripe Chebyshev series approximation: " << chebyshev << std::endl;
  std::cout << "Relative error (%): " << 100 * abs((chebyshev - ana) / ana) << std::endl;
  power = myPowerSeries3dClosed.Evaluate_derivate_x(x, y, z);
  std::cout << "Triple Power series approximation: " << power << std::endl;
  std::cout << "Relative error (%): " << 100 * abs((power - ana) / ana) << std::endl;
  assert(IsClose(ana, chebyshev, 10e-3));
  assert(IsClose(ana, power, 10e-3));

  PrintHeader("Closed segments - y derivative.");
  std::cout << "Point: (x, y, z) = (" << x << ", " << y << ", " << z << ")" << std::endl;
  ana = myFunction3d.Evaluate_derivative_y(x, y, z);
  std::cout << "Analytical: " << ana << std::endl;
  chebyshev = myChebyshevSeries3dClosed.Evaluate_derivative_y(x, y, z);
  std::cout << "Tripe Chebyshev series approximation: " << chebyshev << std::endl;
  std::cout << "Relative error (%): " << 100 * abs((chebyshev - ana) / ana) << std::endl;
  power = myPowerSeries3dClosed.Evaluate_derivate_y(x, y, z);
  std::cout << "Triple Power series approximation: " << power << std::endl;
  std::cout << "Relative error (%): " << 100 * abs((power - ana) / ana) << std::endl;
  assert(IsClose(ana, chebyshev, 10e-3));
  assert(IsClose(ana, power, 10e-3));

  PrintHeader("Closed segments - z derivative.");
  std::cout << "Point: (x, y, z) = (" << x << ", " << y << ", " << z << ")" << std::endl;
  ana = myFunction3d.Evaluate_derivative_z(x, y, z);
  std::cout << "Analytical: " << ana << std::endl;
  chebyshev = myChebyshevSeries3dClosed.Evaluate_derivative_z(x, y, z);
  std::cout << "Tripe Chebyshev series approximation: " << chebyshev << std::endl;
  std::cout << "Relative error (%): " << 100 * abs((chebyshev - ana) / ana) << std::endl;
  power = myPowerSeries3dClosed.Evaluate_derivate_z(x, y, z);
  std::cout << "Triple Power series approximation: " << power << std::endl;
  std::cout << "Relative error (%): " << 100 * abs((power - ana) / ana) << std::endl;
  assert(IsClose(ana, chebyshev, 10e-3));
  assert(IsClose(ana, power, 10e-3));
  std::cout << "\n" << std::endl;

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //                                          Mixed half-open and closed segments
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // aij coefficients.
  xmin = 1.;
  xmax = 30;
  ymin = 1.;
  zmin = 1.;
  order_x = 15;
  order_y = 18;
  order_z = 18;
  x = 10;
  y = 15;
  z = 15;
  auto myChebyshevSeries3dYZOpenedXClosed = ChebyshevSeries3dYZOpenedXClosed<double>(&myFunction3d, xmin, xmax, ymin,
                                                                                     zmin, order_x, order_y, order_z);
  myChebyshevSeries3dYZOpenedXClosed.Compute_aijk();

  // From Chebyshev series to power series.
  bijk = myChebyshevSeries3dYZOpenedXClosed.Compute_bijk();
  auto myPowerSeries3dYZOpenedXClosed = PowerSeries3dYZOpenedXClosed<double>(bijk, xmin, xmax, ymin, zmin);

  // Tests.
  PrintHeader("X closed, YZ open - Function.");
  std::cout << "Point: (x, y, z) = (" << x << ", " << y << ", " << z << ")" << std::endl;
  ana = myFunction3d.Evaluate(x, y, z);
  std::cout << "Analytical: " << ana << std::endl;
  chebyshev = myChebyshevSeries3dYZOpenedXClosed.Evaluate(x, y, z);
  std::cout << "Tripe Chebyshev series approximation: " << chebyshev << std::endl;
  std::cout << "Relative error (%): " << 100 * abs((chebyshev - ana) / ana) << std::endl;
  power = myPowerSeries3dYZOpenedXClosed.Evaluate(x, y, z);
  std::cout << "Triple Power series approximation: " << power << std::endl;
  std::cout << "Relative error (%): " << 100 * abs((power - ana) / ana) << std::endl;
  assert(IsClose(ana, chebyshev, 10e-3));
  assert(IsClose(ana, power, 10e-3));

  PrintHeader("X closed, YZ open - x derivative.");
  std::cout << "Point: (x, y, z) = (" << x << ", " << y << ", " << z << ")" << std::endl;
  ana = myFunction3d.Evaluate_derivative_x(x, y, z);
  std::cout << "Analytical: " << ana << std::endl;
  chebyshev = myChebyshevSeries3dYZOpenedXClosed.Evaluate_derivative_x(x, y, z);
  std::cout << "Tripe Chebyshev series approximation: " << chebyshev << std::endl;
  std::cout << "Relative error (%): " << 100 * abs((chebyshev - ana) / ana) << std::endl;
  power = myPowerSeries3dYZOpenedXClosed.Evaluate_derivate_x(x, y, z);
  std::cout << "Triple Power series approximation: " << power << std::endl;
  std::cout << "Relative error (%): " << 100 * abs((power - ana) / ana) << std::endl;
  assert(IsClose(ana, chebyshev, 10e-3));
  assert(IsClose(ana, power, 10e-3));

  PrintHeader("X closed, YZ open - y derivative.");
  std::cout << "Point: (x, y, z) = (" << x << ", " << y << ", " << z << ")" << std::endl;
  ana = myFunction3d.Evaluate_derivative_y(x, y, z);
  std::cout << "Analytical: " << ana << std::endl;
  chebyshev = myChebyshevSeries3dYZOpenedXClosed.Evaluate_derivative_y(x, y, z);
  std::cout << "Tripe Chebyshev series approximation: " << chebyshev << std::endl;
  std::cout << "Relative error (%): " << 100 * abs((chebyshev - ana) / ana) << std::endl;
  power = myPowerSeries3dYZOpenedXClosed.Evaluate_derivate_y(x, y, z);
  std::cout << "Triple Power series approximation: " << power << std::endl;
  std::cout << "Relative error (%): " << 100 * abs((power - ana) / ana) << std::endl;
  assert(IsClose(ana, chebyshev, 10e-3));
  assert(IsClose(ana, power, 10e-3));

  PrintHeader("X closed, YZ open - z derivative.");
  std::cout << "Point: (x, y, z) = (" << x << ", " << y << ", " << z << ")" << std::endl;
  ana = myFunction3d.Evaluate_derivative_z(x, y, z);
  std::cout << "Analytical: " << ana << std::endl;
  chebyshev = myChebyshevSeries3dYZOpenedXClosed.Evaluate_derivative_z(x, y, z);
  std::cout << "Tripe Chebyshev series approximation: " << chebyshev << std::endl;
  std::cout << "Relative error (%): " << 100 * abs((chebyshev - ana) / ana) << std::endl;
  power = myPowerSeries3dYZOpenedXClosed.Evaluate_derivate_z(x, y, z);
  std::cout << "Triple Power series approximation: " << power << std::endl;
  std::cout << "Relative error (%): " << 100 * abs((power - ana) / ana) << std::endl;
  assert(IsClose(ana, chebyshev, 10e-3));
  assert(IsClose(ana, power, 10e-3));

  std::cout << "\n" << std::endl;

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //                                          Mixed half-open and closed segments
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // aij coefficients.
  xmin = 1.;
  xmax = 30;
  ymin = 1.;
  ymax = 30.;
  zmin = 1.;
  order_x = 15;
  order_y = 15;
  order_z = 18;
  x = 10;
  y = 10;
  z = 15;
  auto myChebyshevSeries3dZOpenedXYClosed = ChebyshevSeries3dZOpenedXYClosed<double>(&myFunction3d, xmin, xmax, ymin,
                                                                                     ymax, zmin, order_x, order_y,
                                                                                     order_z);
  myChebyshevSeries3dZOpenedXYClosed.Compute_aijk();

  // From Chebyshev series to power series.
  bijk = myChebyshevSeries3dZOpenedXYClosed.Compute_bijk();
  auto myPowerSeries3dZOpenedXYClosed = PowerSeries3dZOpenedXYClosed<double>(bijk, xmin, xmax, ymin, ymax, zmin);

  // Tests.
  PrintHeader("XY closed, Z open - Function.");
  std::cout << "Point: (x, y, z) = (" << x << ", " << y << ", " << z << ")" << std::endl;
  ana = myFunction3d.Evaluate(x, y, z);
  std::cout << "Analytical: " << ana << std::endl;
  chebyshev = myChebyshevSeries3dZOpenedXYClosed.Evaluate(x, y, z);
  std::cout << "Tripe Chebyshev series approximation: " << chebyshev << std::endl;
  std::cout << "Relative error (%): " << 100 * abs((chebyshev - ana) / ana) << std::endl;
  power = myPowerSeries3dZOpenedXYClosed.Evaluate(x, y, z);
  std::cout << "Triple Power series approximation: " << power << std::endl;
  std::cout << "Relative error (%): " << 100 * abs((power - ana) / ana) << std::endl;
  assert(IsClose(ana, chebyshev, 10e-3));
  assert(IsClose(ana, power, 10e-3));

  PrintHeader("XY closed, Z open - x derivative.");
  std::cout << "Point: (x, y, z) = (" << x << ", " << y << ", " << z << ")" << std::endl;
  ana = myFunction3d.Evaluate_derivative_x(x, y, z);
  std::cout << "Analytical: " << ana << std::endl;
  chebyshev = myChebyshevSeries3dZOpenedXYClosed.Evaluate_derivative_x(x, y, z);
  std::cout << "Tripe Chebyshev series approximation: " << chebyshev << std::endl;
  std::cout << "Relative error (%): " << 100 * abs((chebyshev - ana) / ana) << std::endl;
  power = myPowerSeries3dZOpenedXYClosed.Evaluate_derivate_x(x, y, z);
  std::cout << "Triple Power series approximation: " << power << std::endl;
  std::cout << "Relative error (%): " << 100 * abs((power - ana) / ana) << std::endl;
  assert(IsClose(ana, chebyshev, 10e-3));
  assert(IsClose(ana, power, 10e-3));

  PrintHeader("XY closed, Z open - y derivative.");
  std::cout << "Point: (x, y, z) = (" << x << ", " << y << ", " << z << ")" << std::endl;
  ana = myFunction3d.Evaluate_derivative_y(x, y, z);
  std::cout << "Analytical: " << ana << std::endl;
  chebyshev = myChebyshevSeries3dZOpenedXYClosed.Evaluate_derivative_y(x, y, z);
  std::cout << "Tripe Chebyshev series approximation: " << chebyshev << std::endl;
  std::cout << "Relative error (%): " << 100 * abs((chebyshev - ana) / ana) << std::endl;
  power = myPowerSeries3dZOpenedXYClosed.Evaluate_derivate_y(x, y, z);
  std::cout << "Triple Power series approximation: " << power << std::endl;
  std::cout << "Relative error (%): " << 100 * abs((power - ana) / ana) << std::endl;
  assert(IsClose(ana, chebyshev, 10e-3));
  assert(IsClose(ana, power, 10e-3));

  PrintHeader("XY closed, Z open - z derivative.");
  std::cout << "Point: (x, y, z) = (" << x << ", " << y << ", " << z << ")" << std::endl;
  ana = myFunction3d.Evaluate_derivative_z(x, y, z);
  std::cout << "Analytical: " << ana << std::endl;
  chebyshev = myChebyshevSeries3dZOpenedXYClosed.Evaluate_derivative_z(x, y, z);
  std::cout << "Tripe Chebyshev series approximation: " << chebyshev << std::endl;
  std::cout << "Relative error (%): " << 100 * abs((chebyshev - ana) / ana) << std::endl;
  power = myPowerSeries3dZOpenedXYClosed.Evaluate_derivate_z(x, y, z);
  std::cout << "Triple Power series approximation: " << power << std::endl;
  std::cout << "Relative error (%): " << 100 * abs((power - ana) / ana) << std::endl;
  assert(IsClose(ana, chebyshev, 10e-3));
  assert(IsClose(ana, power, 10e-3));

}
