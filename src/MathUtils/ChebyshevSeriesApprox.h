//
// Created by pierre-yves on 30/11/2020.
//

#ifndef MATHUTILS_CHEBYSHEVSERIESAPPROX_H
#define MATHUTILS_CHEBYSHEVSERIESAPPROX_H

#include <memory>
#include "Matrix.h"
#include "Functions.h"

namespace mathutils {

/**
* Class for handling the 2d function for the double Chebyshev series approximation.
*/
  template<typename T>
  class Function2d {

   public:

    /// This function evaluates the function at the point (x, y).
    virtual T Evaluate(const double &x, const double &y) const = 0;

  };

/**
 * Class for computing the double Chebyshev series approximation of a function.
 */
  template<typename T>
  class ChebyshevSeriesApprox {

   public:

    /// Contructor of the class.
    ChebyshevSeriesApprox(Function2d<T> *F, const double &xmin, const double &xmax, const double &ymin,
                          const double &ymax, const int &order_x, const int &order_y) : m_function(F), m_x_min(xmin),
                          m_x_max(xmax), m_y_min(ymin), m_y_max(ymax), m_order_x(order_x), m_order_y(order_y) {
      m_aij = MatrixMN<T>::Zero(order_x + 1, order_y + 1);
    }

    /// This method computes the coefficients aij.
    void Computation_aij() {

      // Number of coefficients in x and y.
      int Nx = m_order_x + 1;
      int Ny = m_order_y + 1;

      // x abscissa.
      std::vector<double> x_tilde; // In [-1,1].
      std::vector<double> x; // In [xmin,xmax].

      // x.
      for (int r = 0; r <= m_order_x; ++r) {
        double tmp = cos(MU_PI_2 * (2. * r + 1.) / (m_order_x + 1.));
        x_tilde.push_back(tmp);
        x.push_back(0.5 * (m_x_max - m_x_min) * tmp + 0.5 * (m_x_max + m_x_min));
      }

      // y abscissa.
      std::vector<double> y_tilde; // In [-1,1].
      std::vector<double> y; // In [ymin,ymax].

      // y.
      for (int s = 0; s <= m_order_y; ++s) {
        double tmp = cos(MU_PI_2 * (2. * s + 1.) / (m_order_y + 1.));
        y_tilde.push_back(tmp);
        y.push_back(0.5 * (m_y_max - m_y_min) * tmp + 0.5 * (m_y_max + m_y_min));
      }

      // aij.
      for (int i = 0; i <= m_order_x; ++i) {
        for (int j = 0; j <= m_order_y; ++j) {
          for (int r = 0; r <= m_order_x; ++r) {
            for (int s = 0; s <= m_order_y; ++s) {
              m_aij(i, j) += this->m_function->Evaluate(x.at(r), y.at(s)) * Chebyshev_polynomial(i, x_tilde.at(r))
                             * Chebyshev_polynomial(j, y_tilde.at(s));
              //TODO: Use Chebyshev_polynomial_next instead of Chebyshev_polynomial.
            }
          }
          m_aij(i, j) *= 4. / ((m_order_x + 1.) * (m_order_y + 1.));
          if (i == 0 and j == 0) {
            m_aij(i, j) /= 4.;
          } else if ((i == 0 and j != 0) or (i != 0 and j == 0)) {
            m_aij(i, j) /= 2.;
          }
        }
      }
    }

    /// This method computes the double Chebyshev series approximation.
    T Evaluate(const double &x, const double &y) const {

      double normal_x = (2. / (m_x_max - m_x_min)) * (x - 0.5 * (m_x_max + m_x_min));
      double normal_y = (2. / (m_y_max - m_y_min)) * (y - 0.5 * (m_y_max + m_y_min));

      T result = 0.;
      for (int i = 0; i <= m_order_x; ++i) {
        for (int j = 0; j <= m_order_y; ++j) {
          result += m_aij(i, j) * Chebyshev_polynomial(i, normal_x) * Chebyshev_polynomial(j, normal_y);
        }
      }
      return result;

    }

    /// This method computes the x-derivative double Chebyshev series approximation.
    T Evaluate_derivative_x(const double &x, const double &y) const {

      double normal_x = (2. / (m_x_max - m_x_min)) * (x - 0.5 * (m_x_max + m_x_min));
      double normal_y = (2. / (m_y_max - m_y_min)) * (y - 0.5 * (m_y_max + m_y_min));

      T result = 0.;
      for (int i = 0; i <= m_order_x; ++i) {
        for (int j = 0; j <= m_order_y; ++j) {
          result += m_aij(i, j) * Chebyshev_polynomial_derivative(i, normal_x) * Chebyshev_polynomial(j, normal_y);
        }
      }
      result *= (2. / (m_x_max - m_x_min));
      return result;

    }

    /// This method computes the y-derivative double Chebyshev series approximation.
    T Evaluate_derivative_y(const double &x, const double &y) const {

      double normal_x = (2. / (m_x_max - m_x_min)) * (x - 0.5 * (m_x_max + m_x_min));
      double normal_y = (2. / (m_y_max - m_y_min)) * (y - 0.5 * (m_y_max + m_y_min));

      T result = 0.;
      for (int i = 0; i <= m_order_x; ++i) {
        for (int j = 0; j <= m_order_y; ++j) {
          result += m_aij(i, j) * Chebyshev_polynomial(i, normal_x) * Chebyshev_polynomial_derivative(j, normal_y);
        }
      }
      result *= (2. / (m_y_max - m_y_min));
      return result;

    }

   private:

    /// 2d function to be approximated.
    Function2d<T> *m_function;

    /// Order of the series approximation for x.
    int m_order_x;

    /// Order of the series approximation for y.
    int m_order_y;

    /// x minimum.
    double m_x_min;

    /// x maximum.
    double m_x_max;

    /// y minimum.
    double m_y_min;

    /// y maximum.
    double m_y_max;

    /// aij coefficients.
    MatrixMN <T> m_aij;

  };

}

#endif //MATHUTILS_CHEBYSHEVSERIESAPPROX_H
