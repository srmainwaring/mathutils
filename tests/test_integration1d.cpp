//
// Created by frongere on 18/12/17.
//

#include "MathUtils/MathUtils.h"

using namespace mathutils;

void PrintHeader(std::string title) {
  std::cout << "\n=====================================================================" << std::endl;
  std::cout << "    " << title << std::endl;
  std::cout << "=====================================================================" << std::endl;
}

// Definition of a function to integrate.
double myFunction(double x) {
  if (x == 0.) {
      return 1.;
  } else {
      return sin(x) / x;
  }
}

int main(int argc, char* argv[]) {

  // This test checks the numerical integration of a one-variable function.

  PrintHeader("Computation of the integral of sinc(x) between 0 and 1.");
  std::cout << "Analytical: 0.946083" << std::endl;

  // Defining the integration.
  // Here, by default we use the trapezoidal method.
  auto myIntegrator = Integrate1d<double>(myFunction, 0, 1, 1000);
  std::cout << "Integration by a function: " << myIntegrator.Get() << std::endl;

  // Defining a vector.
  auto x = linspace<double>(0, 1, 1000);
  std::vector<double> y;
  for (auto val : x) {
      y.push_back(myFunction(val));
  }
  auto myIntegratorVect = Integrate1d<double>(y, 0, 1, 1000);
  std::cout << "Integration by a vector: " << myIntegratorVect.Get() << std::endl;

  return 0;
}