//
// Created by pierre-yves on 31/10/20.
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

  // This test checks the computation of the special functions.

  // ========================================================================
  //                          Bessel function
  // ========================================================================
  PrintHeader("Bessel function");
  std::cout << "" << std::endl;

  // Source: https://en.cppreference.com/w/cpp/numeric/special_functions/cyl_bessel_j
  PrintInfo("Comparison of the value returned by the method and an exact value for n = 0 and x = 1.2345");
  double value_J0 = Cyl_Bessel_first_kind<int, double, double>(0, 1.2345);
  std::cout << "Method Cyl_Bessel_first_kind: " << value_J0 << std::endl;
  std::cout << "Exact value: 0.653792" << std::endl;
  assert(IsClose(value_J0, 0.653792));

  // ========================================================================
  //                       Legendre polynomials
  // ========================================================================
  PrintHeader("Legendre polynomials");
  std::cout << "" << std::endl;

  // Source: https://en.cppreference.com/w/cpp/numeric/special_functions/legendre
  PrintInfo("Comparison of the value returned by the method and an exact value for n = 3 and x = 0.25");
  double value_Legendre_P3 = Legendre_polynomial(3, 0.25);
  std::cout << "Method Legendre: " << value_Legendre_P3 << std::endl;
  std::cout << "Exact value: " << 0.5*(5*0.25*0.25*0.25 - 3*0.25) << std::endl; // P_3(x) = 1/2*(5x^3 - 3x)
  assert(IsClose(value_Legendre_P3, 0.5*(5*0.25*0.25*0.25 - 3*0.25)));

  // ========================================================================
  //                 Derivative of the Legendre polynomials
  // ========================================================================
  PrintHeader("Derivative of the Legendre polynomials");
  std::cout << "" << std::endl;

  PrintInfo("Comparison of the value returned by the method and an exact value for n = 3 and x = 0.25");
  double value_derivative_Legendre_P3 = Legendre_polynomial_derivative(3, 0.25);
  std::cout << "Method Legendre_derivative: " << value_derivative_Legendre_P3 << std::endl;
  std::cout << "Exact value: " << 0.5*(15*0.25*0.25 - 3) << std::endl; // 1/2*(15x^2 - 3)
  assert(IsClose(value_derivative_Legendre_P3, 0.5*(15*0.25*0.25 - 3)));

  // ========================================================================
  //                            Factorial
  // ========================================================================
  PrintHeader("Factorial");
  std::cout << "" << std::endl;

  PrintInfo("Comparison of the value returned by the method and an exact value for n = 9");
  double value_factorial = Factorial(9);
  std::cout << "Method Factorial: " << value_factorial << std::endl;
  std::cout << "Exact value: " << 362880. << std::endl; // 9! = 362880
  assert(IsClose(value_factorial, 362880.));

}