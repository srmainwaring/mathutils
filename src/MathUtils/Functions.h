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

  // This function returns the exponential integral Ei at the point x.
  template <class T>
  T Ei(const T x){
    return boost::math::expint(x);
  }

  // This function returns the exponential integral Ei at the point x evaluated by Chebyshev and continued fraction approximations.
  template <class T>
  T Ei_approximation(const T x){

    // This approximation comes from the following article.
    // Cody W. J. and Thacher H. C. Chebyshev Approximations for the Exponential Integral Ei(x), Mathematics of Computation, 23(106): 289-303, 1969.
    // doi: 10.2307/2004423.

    T Ei = 0.;

    if(IsClose(x, 0.)) {
      std::cout << "Ei_Chebyshev_approximation: Ei is not defined for x = 0." << std::endl;
      exit(0);
    }

    if(x > 0. and x < 6.){

      // Parameters.
      int n = 9;
      double x0 = 0.372507410781366634461991866580;
      std::vector<double> pj = {-4.1658081333604994241879e11,
                                1.2177698136199594677580e10,
                                -2.5301823984599019348858e10,
                                3.1984354235237738511048e8,
                                -3.5377809694431133484800e8,
                                -3.1398660864247265862050e5,
                                -1.4299841572091610380064e6,
                                -1.4287072500197005777376e4,
                                -1.2831220659262000678155e3,
                                -1.2963702602474830028590e1};
      std::vector<double> qj = {-1.7934749837151009723371e11,
                                9.8900934262481749439886e10,
                                -2.8986272696554495342658e10,
                                5.4229617984472955011862e9,
                                -7.0108568774215954065376e8,
                                6.4698830956576428587653e7,
                                -4.2648434812177161405483e6,
                                1.9418469440759880361415e5,
                                -5.5648470543369082846819e3,
                                7.688671875000000000000e1};
      assert (pj.size() == qj.size());

      // Evaluation.
      Ei += log(x / x0);
      T numerator = 0.5 * pj.at(0) * Chebyshev_polynomial_shifted<int, T>(0, x / 6.);
      T denominator = 0.5 * qj.at(0) * Chebyshev_polynomial_shifted<int, T>(0, x / 6.);
      for (unsigned int j = 1; j <= n; ++j) {
        numerator += pj.at(j) * Chebyshev_polynomial_shifted<int, T>(j, x / 6.);
        denominator += qj.at(j) * Chebyshev_polynomial_shifted<int, T>(j, x / 6.);
      }
      Ei += (x - x0) * (numerator / denominator);

    }
    else if(x >= 6. and x < 12.) {

      // Parameters.
      int n = 9;
      std::vector<double> alphaj = {9.981193787537396413219e-1,
                                    9.565134591978630774217e0,
                                    -3.988850730390541057912e0,
                                    1.120011024227297451523e1,
                                    -3.015761863840593359165e1,
                                    1.945603779539281810439e1,
                                    1.052976392459015155422e1,
                                    -2.421106956980653511550e1,
                                    -2.378372882815725244124e0,
                                    -2.645977793077147237806e0};
      std::vector<double> betaj = {1.249884822712447891440e0,
                                   -2.369210235636181001661e2,
                                   4.731097187816050252967e2,
                                   2.852397548119248700147e1,
                                   7.608194509086645763123e2,
                                   -8.791401054875438925029e0,
                                   3.697412299772985940785e2,
                                   4.644185932583286942650e0,
                                   1.598517957704779356479e-4};
      assert (alphaj.size() == (betaj.size() + 1));

      // Evaluation.
      T fraction = 0.;
      for (unsigned int j = n; j >= 1; --j) {
        fraction = (betaj.at(j - 1) / (alphaj.at(j) + x + fraction));
      }
      fraction += alphaj.at(0); // First term.
      Ei = (exp(x) / x) * fraction;

    }
    else if(x >= 12. and x < 24.) {

      // Parameters.
      int n = 9;
      std::vector<double> alphaj = {9.999933106160568739091e-1,
                                    -1.845086232391278674524e0,
                                    2.652575818452799819855e1,
                                    2.495487730402059440626e1,
                                    -3.323612579343962284333e1,
                                    -9.134835699998742552432e-1,
                                    -2.105740799548040450394e1,
                                    -1.000641913989284829961e1,
                                    -1.860092121726437582253e1,
                                    -1.647721172463463140042e0};
      std::vector<double> betaj = {1.001533852045342697818e0,
                                   -1.093556195391091243924e1,
                                   1.991004470817742470726e2,
                                   1.192832423968601006985e3,
                                   4.429413178337928401161e1,
                                   2.538819315630708031713e2,
                                   5.994932325667407355255e1,
                                   6.403800405352415551324e1,
                                   9.792403599217290296840e1};
      assert (alphaj.size() == (betaj.size() + 1));

      // Evaluation.
      T fraction = 0.;
      for (unsigned int j = n; j >= 1; --j) {
        fraction = (betaj.at(j - 1) / (alphaj.at(j) + x + fraction));
      }
      fraction += alphaj.at(0); // First term.
      Ei = (exp(x) / x) * fraction;

    }
    else { // x >= 24.

      // Parameters.
      int n = 9;
      std::vector<double> alphaj = {1.00000000000000485503e0,
                                    -3.00000000320981265753e0,
                                    -5.00006640413131002475e0,
                                    -7.06810977895029358836e0,
                                    -1.52856623636929636839e1,
                                    -7.63147701620253630855e0,
                                    -2.79798528624305389340e1,
                                    -1.81949664929868906455e1,
                                    -2.23127670777632409550e2,
                                    1.75338801265465972390e2};
      std::vector<double> betaj = {1.99999999999048104167e0,
                                   -2.99999894040324959612e0,
                                   -7.99243595776339741065e0,
                                   -1.20187763547154743238e1,
                                   7.04831847180424675988e1,
                                   1.17179220502086455287e2,
                                   1.37790390235747998793e2,
                                   3.97277109100414518365e0,
                                   3.97845977167414720840e4};
      assert (alphaj.size() == (betaj.size() + 1));

      // Evaluation.
      T fraction = 0.;
      for (unsigned int j = n; j >= 1; --j) {
        fraction = (betaj.at(j - 1) / (alphaj.at(j) + x + fraction));
      }
      fraction += alphaj.at(0); // First term.
      Ei = (exp(x) / x) * (1. + (1. / x) * fraction);

    }

    return Ei;

  }

  // This function returns the quantity exp(-x)*Ei(x).
  template <class T>
  T expEi(const T x){

    T result = 0.;

    if(x < 24.){
      result = exp(-x) * Ei(x);
    }
    else{ // x >= 24.

      // Parameters.
      int n = 9;
      std::vector<double> alphaj = {1.00000000000000485503e0,
                                    -3.00000000320981265753e0,
                                    -5.00006640413131002475e0,
                                    -7.06810977895029358836e0,
                                    -1.52856623636929636839e1,
                                    -7.63147701620253630855e0,
                                    -2.79798528624305389340e1,
                                    -1.81949664929868906455e1,
                                    -2.23127670777632409550e2,
                                    1.75338801265465972390e2};
      std::vector<double> betaj = {1.99999999999048104167e0,
                                   -2.99999894040324959612e0,
                                   -7.99243595776339741065e0,
                                   -1.20187763547154743238e1,
                                   7.04831847180424675988e1,
                                   1.17179220502086455287e2,
                                   1.37790390235747998793e2,
                                   3.97277109100414518365e0,
                                   3.97845977167414720840e4};
      assert (alphaj.size() == (betaj.size() + 1));

      // Evaluation.
      T fraction = 0.;
      for (unsigned int j = n; j >= 1; --j) {
        fraction = (betaj.at(j - 1) / (alphaj.at(j) + x + fraction));
      }
      fraction += alphaj.at(0); // First term.
      result = (1. / x) * (1. + (1. / x) * fraction);
    }

    return result;

  }

}  // end namespace mathutils

#endif //MATHUTILS_FUNCTIONS_H
