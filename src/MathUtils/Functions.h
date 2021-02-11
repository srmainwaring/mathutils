//
// Created by frongere on 28/12/17.
//

#ifndef MATHUTILS_FUNCTIONS_H
#define MATHUTILS_FUNCTIONS_H

#include "StdVector.h"
#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/special_functions/legendre.hpp>
#include <boost/math/special_functions/factorials.hpp>
#include <boost/math/special_functions/expint.hpp>
#include <boost/math/special_functions/chebyshev.hpp>
#include "VectorGeneration.h"

namespace mathutils {

  template <class Real>
  std::vector<Real> Sinc(const std::vector<Real>& x) {
      unsigned long n = x.size();
      std::vector<Real> out(n);
      for (unsigned long i=0; i<n; i++) {
          if (x[i] == 0.) {
              out[i] = 1.;
          } else {
              out[i] = sin(x[i]) / x[i];
          }
      }
      return out;
  }

  template <class Real>
  std::vector<Real> Sinc(const Real xmin, const Real xmax, const unsigned int n) {
      auto x = linspace<Real>(xmin, xmax, n);
      return Sinc(x);
  }

  template <class Real>
  std::vector<Real> Heaviside(const std::vector<Real>& x, const Real x_step) {
      unsigned long n = x.size();

      std::vector<Real> out(n);
      Real xval;
      for (unsigned long i=0; i<n; i++) {
          xval = x[i];
          if (xval < x_step) {
              out[i] = 0.;
          } else if (xval > x_step) {
              out[i] = 1.;
          } else {
              out[i] = 0.5;
          }
      }
      return out;
  }

  template <class Real>
  std::vector<Real> Door(const std::vector<Real>& x, const Real x_min, const Real x_max, const Real height=1.) {
      auto out = Heaviside(x, x_min);
      out = Minus(out, Heaviside<Real>(x, x_max));
      Scale<Real>(out, height);
      return out;
  }

  template <class Real>
  std::vector<Real> Dirac(const std::vector<Real>& x, const Real x_impulse, const Real height=1.) {
      unsigned long n = x.size();

      Real dx_2 = 0.5 * (x[1] - x[0]);

      std::vector<Real> out(n);
      for (unsigned long i=0; i<n; i++) {
          if (abs(x[i] - x_impulse) < dx_2) { // FIXME: le critere n'est pas robuste
              out[i] = height;
          } else {
              out[i] = 0.;
          }
      }
      return out;
  }

  template <class Real>
  std::vector<Real> Ramp(const std::vector<Real>& x,
                         const Real x_min, const Real x_max,
                         const Real y_left, const Real y_right) {

      unsigned long n = x.size();

      assert(x_min <= x_max);
      assert(x_min >= x[0] && x_max <= x[n-1]);

      Real a = (y_right - y_left) / (x_max - x_min);
      Real b = y_left - a*x_min;

      std::vector<Real> out(n);
      Real xi;
      for (unsigned long i=0; i<n; i++) {
          xi = x[i];
          if (xi <= x_min) {
              out[i] = y_left;
          } else if (x[i] <= x_max) {
              out[i] = a*xi + b;
          } else {
              out[i] = y_right;
          }
      }

      return out;
  }

  // This function returns the Bessel function of the first kind of order "order" at the point x.
  template <class T1, class T2, class Tresults>
  Tresults Cyl_Bessel_first_kind(const T1& order, const T2 x){
      return boost::math::cyl_bessel_j(order, x);
  }

  // This function returns the Bessel function of the second kind of order "order" at the point x.
  template <class T1, class T2, class Tresults>
  Tresults Cyl_Bessel_second_kind(const T1& order, const T2 x){
    return boost::math::cyl_neumann(order, x);
  }

  // This function returns the modified Bessel function of the first kind of order "order" at the point x.
  template <class T1, class T2, class Tresults>
  Tresults Cyl_modified_Bessel_first_kind(const T1& order, const T2 x){
    return boost::math::cyl_bessel_i(order, x);
  }

  // This function returns the modified Bessel function of the second kind of order "order" at the point x.
  template <class T1, class T2, class Tresults>
  Tresults Cyl_modified_Bessel_second_kind(const T1& order, const T2 x){
    return boost::math::cyl_bessel_k(order, x);
  }

  // This function returns the Legendre polynomial of order "order" at the point x.
  template <class Integer, class Real>
  Real Legendre_polynomial(const Integer& order, const Real x){
      return boost::math::legendre_p(order, x);
  }

  // This function returns the derivative of the Legendre polynomial of order "order" at the point x.
  template <class Integer, class Real>
  Real Legendre_polynomial_derivative(const Integer& order, const Real x){
    return boost::math::legendre_p_prime(order, x);
  }

  // This function returns the factorial of the integer n.
  template <class Integer, class Real>
  Real Factorial(const Integer n){
    return boost::math::factorial<double>(n);
  }
    // This function returns the exponential integral Ei at the point x.
    template <class T>
    T Ei(const T x){
      return boost::math::expint(x);
    }

  // This function returns the zero-order Struve function at the point x.
  template <class T>
  T Struve_zero_order(const T x){
    T Struve = 0.;
    if(x >= 0. and x <= 3.) {
      T arg = x / 3.;
      std::vector<double> coef_a = {1.909859164, -1.909855001, 0.687514637, -0.126164557, 0.013828813, -0.000876918};
      for (unsigned int j = 1; j <= 6; ++j) {
        Struve += coef_a.at(j - 1) * pow(arg, 2 * j - 1);
      }
    }
    else if(x > 3.) {
      T arg = 3. / x;
      std::vector<double> coef_b = {0.99999906, 4.77228920, 3.85542044, 0.32303607};
      std::vector<double> coef_c = {1., 4.88331068, 4.28957333, 0.52120508};
      T numerator = 0.;
      T denominator = 0.;
      for (unsigned int j = 0; j <= 3; ++j) {
        numerator += coef_b.at(j) * pow(arg, 2 * j);
        denominator += coef_c.at(j) * pow(arg, 2 * j);
      }
      Struve = Cyl_Bessel_second_kind<int, T, T>(0, x) + (2. * numerator / (MU_PI * x * denominator));
    }
    else {
      std::cout << "Zero-order Struve function is only defined for a positive argument." << std::endl;
      exit(0);
    }

    return Struve;
  }

  // This function returns the first-order Struve function at the point x.
  template <class T>
  T Struve_first_order(const T x){
    T Struve = 0.;
    if(x >= 0. and x <= 3.) {
      T arg = x / 3.;
      std::vector<double> coef_d = {1.909859286, -1.145914713, 0.294656958, -0.042070508, 0.003785727, -0.000207183};
      for (unsigned int j = 1; j <= 6; ++j) {
        Struve += coef_d.at(j - 1) * pow(arg, 2 * j);
      }
    }
    else if(x > 3.) {
      T arg = 3. / x;
      std::vector<double> coef_e = {1.00000004, 3.92205313, 2.64893033, 0.27450895};
      std::vector<double> coef_f = {1., 3.81095112, 2.26216956, 0.10885141};
      T numerator = 0.;
      T denominator = 0.;
      for (unsigned int j = 0; j <= 3; ++j) {
        numerator += coef_e.at(j) * pow(arg, 2 * j);
        denominator += coef_f.at(j) * pow(arg, 2 * j);
      }
      Struve = Cyl_Bessel_second_kind<int, T, T>(1, x) + (2. * numerator / (MU_PI * denominator));
    }
    else {
      std::cout << "First-order Struve function is only defined for a positive argument." << std::endl;
      exit(0);
    }

    return Struve;
  }

  // This function returns the derivative of the zero-order Struve function at the point x.
  template <class T>
  T Struve_zero_order_derivative(const T x){
    return (MU_2_PI - Struve_first_order<T>(x));
  }

  // This function returns the Chebyshev polynomial of order "order" at the point x.
  template <class Integer, class Real>
  Real Chebyshev_polynomial(const Integer& order, const Real x){
    return boost::math::chebyshev_t(order, x);
  }

  // This function returns a shifted Chebyshev polynomial of order "order" at the point x as used in
  // Cody and Tchacher 1969 in the function Ei_Chebyshev_approximation.
  template <class Integer, class Real>
  Real Chebyshev_polynomial_shifted(const Integer& order, const Real x){
    // Ti_star(x) = Ti(2x - 1).
    return boost::math::chebyshev_t(order, 2. * x - 1.);
  }

  // This function returns the next Chebyshev polynomial at the point x, given the two previous polynomials.
  template <class Real>
  Real Chebyshev_polynomial_next(const Real x, const Real Tn, const Real Tn_1){
    return boost::math::chebyshev_next(x, Tn, Tn_1);
  }

  // This function returns the derivative of the Chebyshev polynomial of order "order" at the point x.
  template <class Integer, class Real>
  Real Chebyshev_polynomial_derivative(const Integer& order, const Real x){
    return boost::math::chebyshev_t_prime(order, x);
  }

}  // end namespace mathutils

#endif //MATHUTILS_FUNCTIONS_H
