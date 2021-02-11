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
  //                          Bessel functions
  // ========================================================================
  PrintHeader("Bessel function - First kind");
  std::cout << "" << std::endl;

  // Source: https://en.cppreference.com/w/cpp/numeric/special_functions/cyl_bessel_j
  PrintInfo("Comparison of the value returned by the method and an exact value for n = 0 and x = 1.2345");
  double value_J0 = Cyl_Bessel_first_kind<int, double, double>(0, 1.2345);
  std::cout << "Method Cyl_Bessel_first_kind: " << value_J0 << std::endl;
  std::cout << "Exact value: 0.653792" << std::endl;
  assert(IsClose(value_J0, 0.653792));

  PrintHeader("Bessel function - Second kind");
  std::cout << "" << std::endl;

  // Source: https://en.cppreference.com/w/cpp/numeric/special_functions/cyl_neumann
  PrintInfo("Comparison of the value returned by the method and an exact value for n = 0.5 and x = 1.2345");
  double value_Y05 = Cyl_Bessel_second_kind<double, double, double>(0.5, 1.2345);
  std::cout << "Method Cyl_Bessel_second_kind: " << value_Y05 << std::endl;
  std::cout << "Exact value: -0.236973" << std::endl;
  assert(IsClose(value_Y05, -0.236973));

  PrintHeader("Modified Bessel function - First kind");

  // Source: https://en.cppreference.com/w/cpp/numeric/special_functions/cyl_bessel_i
  PrintInfo("Comparison of the value returned by the method and an exact value for n = 0 and x = 1.2345");
  double value_I0 = Cyl_modified_Bessel_first_kind<double, double, double>(0, 1.2345);
  std::cout << "Method Cyl_modified_Bessel_first_kind: " << value_I0 << std::endl;
  std::cout << "Exact value: 1.41886" << std::endl;
  assert(IsClose(value_I0, 1.41886));

  PrintHeader("Modified Bessel function - Second kind");

  // Source: https://en.cppreference.com/w/cpp/numeric/special_functions/cyl_bessel_k
  PrintInfo("Comparison of the value returned by the method and an exact value for n = 0.5 and x = 1.2345");
  double value_K05 = Cyl_modified_Bessel_second_kind<double, double, double>(0.5, 1.2345);
  std::cout << "Method Cyl_modified_Bessel_second_kind: " << value_K05 << std::endl;
  std::cout << "Exact value: 0.32823" << std::endl;
  assert(IsClose(value_K05, 0.32823));

  // ========================================================================
  //                       Legendre polynomials
  // ========================================================================
  PrintHeader("Legendre polynomials");
  std::cout << "" << std::endl;

  // Source: https://en.cppreference.com/w/cpp/numeric/special_functions/legendre
  PrintInfo("Comparison of the value returned by the method and an exact value for n = 3 and x = 0.25");
  double value_Legendre_P3 = Legendre_polynomial<int, double>(3, 0.25);
  std::cout << "Method Legendre_polynomial: " << value_Legendre_P3 << std::endl;
  std::cout << "Exact value: " << 0.5*(5*0.25*0.25*0.25 - 3*0.25) << std::endl; // P_3(x) = 1/2*(5x^3 - 3x)
  assert(IsClose(value_Legendre_P3, 0.5*(5*0.25*0.25*0.25 - 3*0.25)));

  // ========================================================================
  //                 Derivative of the Legendre polynomials
  // ========================================================================
  PrintHeader("Derivative of the Legendre polynomials");
  std::cout << "" << std::endl;

  PrintInfo("Comparison of the value returned by the method and an exact value for n = 3 and x = 0.25");
  double value_derivative_Legendre_P3 = Legendre_polynomial_derivative<int, double>(3, 0.25);
  std::cout << "Method Legendre_derivative: " << value_derivative_Legendre_P3 << std::endl;
  std::cout << "Exact value: " << 0.5*(15*0.25*0.25 - 3) << std::endl; // 1/2*(15x^2 - 3)
  assert(IsClose(value_derivative_Legendre_P3, 0.5*(15*0.25*0.25 - 3)));

  // ========================================================================
  //                            Factorial
  // ========================================================================
  PrintHeader("Factorial");
  std::cout << "" << std::endl;

  PrintInfo("Comparison of the value returned by the method and an exact value for n = 9");
  double value_factorial = Factorial<int, double>(9);
  std::cout << "Method Factorial: " << value_factorial << std::endl;
  std::cout << "Exact value: " << 362880. << std::endl; // 9! = 362880
  assert(IsClose(value_factorial, 362880.));

  // ========================================================================
  //                        Exponential integral
  // ========================================================================
  PrintHeader("Exponential integral");
  std::cout << "" << std::endl;

  // Source: https://en.cppreference.com/w/cpp/numeric/special_functions/expint
  PrintInfo("Comparison of the value returned by the direct method and an exact value for x = 1");
  double value_Ei = Ei<double>(1.);
  std::cout << "Method Ei: " << value_Ei << std::endl;
  std::cout << "Exact value: 1.89512" << std::endl;
  assert(IsClose(value_Ei, 1.89512));

  PrintInfo("Comparison of the value returned by the Chebyshev approximation method and an exact value for x = 1");
  double value_Ei_approx = Ei_approximation<double>(1.);
  std::cout << "Method Ei_Chebyshev_approximation: " << value_Ei_approx << std::endl;
  std::cout << "Exact value: 1.89512" << std::endl;
  assert(IsClose(value_Ei_approx, 1.89512));
  
  PrintInfo("\nComparison of the value returned by the two methods for evaluating Ei");
  PrintInfo("x = 3");
  value_Ei = Ei<double>(3.);
  value_Ei_approx = Ei_approximation<double>(3.);
  std::cout << "Method Ei: " << value_Ei << std::endl;
  std::cout << "Method Ei_Chebyshev_approximation: " << value_Ei_approx << std::endl;
  std::cout << "|Delta| = " << abs(value_Ei - value_Ei_approx) << std::endl;
  assert(IsClose(abs(value_Ei - value_Ei_approx), 0.));
  PrintInfo("\nx = 9");
  value_Ei = Ei<double>(9.);
  value_Ei_approx = Ei_approximation<double>(9.);
  std::cout << "Method Ei: " << value_Ei << std::endl;
  std::cout << "Method Ei_Chebyshev_approximation: " << value_Ei_approx << std::endl;
  std::cout << "|Delta| = " << abs(value_Ei - value_Ei_approx) << std::endl;
  assert(IsClose(abs(value_Ei - value_Ei_approx), 0.));
  PrintInfo("\nx = 18");
  value_Ei = Ei<double>(18.);
  value_Ei_approx = Ei_approximation<double>(18.);
  std::cout << "Method Ei: " << value_Ei << std::endl;
  std::cout << "Method Ei_Chebyshev_approximation: " << value_Ei_approx << std::endl;
  std::cout << "|Delta| = " << abs(value_Ei - value_Ei_approx) << std::endl;
  assert(IsClose(abs(value_Ei - value_Ei_approx), 0.));
  PrintInfo("\nx = 30");
  value_Ei = Ei<double>(30.);
  value_Ei_approx = Ei_approximation<double>(30.);
  std::cout << "Method Ei: " << value_Ei << std::endl;
  std::cout << "Method Ei_Chebyshev_approximation: " << value_Ei_approx << std::endl;
  std::cout << "|Delta| = " << abs(value_Ei - value_Ei_approx) << std::endl;
  assert(IsClose(abs(value_Ei - value_Ei_approx), 0.));

  // ========================================================================
  //                       Chebyshev polynomials
  // ========================================================================
  PrintHeader("Chebyshev polynomials");
  std::cout << "" << std::endl;

  // Source: https://fr.wikipedia.org/wiki/Polyn%C3%B4me_de_Tchebychev
  PrintInfo("Comparison of the value returned by the method and an exact value for n = 5 and x = 0.25");
  double value_Chebyshev_T5_1 = Chebyshev_polynomial<int, double>(5, 0.25);
  std::cout << "Method Chebyshev_polynomial: " << value_Chebyshev_T5_1 << std::endl;
  double value_Chebyshev_T3 = Chebyshev_polynomial<int, double>(3, 0.25);
  double value_Chebyshev_T4 = Chebyshev_polynomial<int, double>(4, 0.25);
  double value_Chebyshev_T5_2 = Chebyshev_polynomial_next<double>(0.25, value_Chebyshev_T4, value_Chebyshev_T3);
  std::cout << "Method Chebyshev_polynomial_next: " << value_Chebyshev_T5_2 << std::endl;
  std::cout << "Exact value: " << 16.*pow(0.25, 5.) - 20.*pow(0.25, 3.) + 5.*0.25 << std::endl; // T_5(x) = 16x^5 - 20x^3 + 5x.
  assert(IsClose(value_Chebyshev_T5_1, 16.*pow(0.25, 5.) - 20.*pow(0.25, 3.) + 5.*0.25));
  assert(IsClose(value_Chebyshev_T5_2, 16.*pow(0.25, 5.) - 20.*pow(0.25, 3.) + 5.*0.25));

  // ========================================================================
  //                   Derivate of Chebyshev polynomials
  // ========================================================================
  PrintHeader("Derivate of Chebyshev polynomials");
  std::cout << "" << std::endl;

  PrintInfo("Comparison of the value returned by the method and an exact value for n = 5 and x = 0.25");
  double value_Chebyshev_derivative_T5 = Chebyshev_polynomial_derivative<int, double>(5, 0.25);
  std::cout << "Method Chebyshev_polynomial_derivative: " << value_Chebyshev_derivative_T5 << std::endl;
  std::cout << "Exact value: " << 80.*pow(0.25, 4.) - 60.*pow(0.25, 2.) + 5 << std::endl; // T'_5(x) = 80x^4 - 60x^2 + 5.
  assert(IsClose(value_Chebyshev_derivative_T5, 80.*pow(0.25, 4.) - 60.*pow(0.25, 2.) + 5));

}