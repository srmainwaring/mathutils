//
// Created by pierre-yves on 03/03/2021.
//

#ifndef MATHUTILS_CHEBYSHEVTOPOWERSERIES2D_H
#define MATHUTILS_CHEBYSHEVTOPOWERSERIES2D_H

#include "Matrix.h"
#include "Polynomial.h"
#include "ChebyshevSeries2d.h"
#include "Functions.h"

namespace mathutils {

/**
* Class for handling and computing the double power series approximation of a function from its Chebyshev series approximation.
*/
  template<typename T>
  class ChebyshevToPowerSeries2dBase {

   public:

    /// Contructor of the class.
    ChebyshevToPowerSeries2dBase(const std::shared_ptr<ChebyshevSeries2dBase<T>> &chebyshev) {
      m_x_min = chebyshev->x_min();
      m_y_min = chebyshev->y_min();
      m_order_x = chebyshev->order_x();
      m_order_y = chebyshev->order_y();

      if(m_order_x > 18 or m_order_y > 18) {
        std::cout << "The conversion Chebyshev to power series have large numerical inaccuracies for important order." << std::endl;
        std::cout << "The maximum order in x and y is 18." << std::endl;
        std::cout << "Order in x: " << m_order_x << std::endl;
        std::cout << "Order in y: " << m_order_y << std::endl;
        exit(0);
      }

      m_pij = MatrixMN<T>::Zero(m_order_x + 1, m_order_y + 1);
      Compute_pij(chebyshev->aij());
    }

    /// This method computes the transformation coefficient.
    double TransformationCoefficient(const int& lj, const int& kj) {
      double lambda;
      int ljkj_plus = lj + kj;
      if (lj == 0 and kj == 0) {
        lambda = 1.;
      } else if (ljkj_plus % 2 == 1) { // (lj + kj) impair.
        lambda = 0.;
      } else { // (lj + kj) pair.
        int floor_plus = floor(0.5 * ljkj_plus);
        int floor_minus = floor(0.5 * (lj - kj));
        lambda = pow(-1, floor_minus) * pow(2., kj - 1.) * lj * (Factorial<int, double>(floor_plus - 1) /
                    (Factorial<int, double>(floor_minus) * Factorial<int, double>(kj)));
      }

      return lambda;

    }

    /// This method computes the pij of the power series from the aij of the Chebyshev series.
    void Compute_pij(const MatrixMN<T> &aij) {

      for (int k1 = 0; k1 <= m_order_x; ++k1) {
        for (int k2 = 0; k2 <= m_order_y; ++k2) {
          for (int l1 = k1; l1 <= m_order_x; ++l1) {
            double lamda_l1k1 = TransformationCoefficient(l1, k1);
            double tmp = 0.;
            for (int l2 = k2; l2 <= m_order_y; ++l2) {
              double lamda_l2k2 = TransformationCoefficient(l2, k2);
              tmp += aij(l1, l2) * lamda_l2k2;
            }
            m_pij(k1, k2) += tmp * lamda_l1k1;
          }
        }
      }

    }

    /// This method computes the double power series approximation.
    T Evaluate(const double &x, const double &y) const {

      // Parameters.
      double xunit = AffineTransformationSegmentToUnit_x(x);
      double yunit = AffineTransformationSegmentToUnit_y(y);

      // Partial sum.
      std::vector<double> qi;
      qi.reserve(m_order_x + 1);
      for (int i = 0; i <= m_order_x; ++i) {
        Eigen::VectorXd tmp = m_pij.row(i);
        std::vector<double> pj(tmp.data(), tmp.data() + tmp.size());
        qi.push_back(Horner<double>(pj, yunit));
      }
      T result = Horner<double>(qi, xunit);

      return result;

    }

    /// This method computes the x-derivative double power series approximation.
    T Evaluate_derivate_x(const double &x, const double &y) const {

      // Parameters.
      double xunit = AffineTransformationSegmentToUnit_x(x);
      double yunit = AffineTransformationSegmentToUnit_y(y);

      // Partial sum.
      std::vector<double> qi;
      qi.reserve(m_order_x + 1);
      for (int i = 0; i <= m_order_x; ++i) {
        Eigen::VectorXd tmp = m_pij.row(i);
        std::vector<double> pj(tmp.data(), tmp.data() + tmp.size());
        qi.push_back(Horner<double>(pj, yunit));
      }
      T result = CoefficientDerivative_x(x) * Horner_derivative<double>(qi, xunit);

      return result;

    }

    /// This method computes the y-derivative double power series approximation.
    T Evaluate_derivate_y(const double &x, const double &y) const {

      // Parameters.
      double xunit = AffineTransformationSegmentToUnit_x(x);
      double yunit = AffineTransformationSegmentToUnit_y(y);

      // Partial sum.
      std::vector<double> qi;
      qi.reserve(m_order_x + 1);
      for (int i = 0; i <= m_order_x; ++i) {
        Eigen::VectorXd tmp = m_pij.row(i);
        std::vector<double> pj(tmp.data(), tmp.data() + tmp.size());
        qi.push_back(Horner_derivative<double>(pj, yunit));
      }
      T result = CoefficientDerivative_y(x) * Horner<double>(qi, xunit);

      return result;

    }

   protected:

    /// This method applied an affine transformation from [-1, 1] to the domain of the approximation for x.
    virtual double AffineTransformationUnitToSegment_x(const double &xunit) const = 0;

    /// This method applied an affine transformation from [-1, 1] to the domain of the approximation for y.
    virtual double AffineTransformationUnitToSegment_y(const double &yunit) const = 0;

    /// This method applied an affine transformation from the domain of the approximation to [-1, 1] for x.
    virtual double AffineTransformationSegmentToUnit_x(const double &xdomain) const = 0;

    /// This method applied an affine transformation from the domain of the approximation to [-1, 1] for y.
    virtual double AffineTransformationSegmentToUnit_y(const double &ydomain) const = 0;

    /// This method gives the coefficient in front of the derivate of the  double Chebychev series wrt x.
    virtual double CoefficientDerivative_x(const double &xdomain) const = 0;

    /// This method gives the coefficient in front of the derivate of the  double Chebychev series wrt y.
    virtual double CoefficientDerivative_y(const double &ydomain) const = 0;

   protected:

    /// Order of the series approximation for x.
    int m_order_x;

    /// Order of the series approximation for y.
    int m_order_y;

    /// x minimum.
    double m_x_min;

    /// y minimum.
    double m_y_min;

    /// Pij coefficients.
    MatrixMN<T> m_pij;

  };

  /**
  * Class for computing the double power series approximation of a function from its Chebyshev series approximation over
   * [xmin, xmax] x [ymin, ymax].
  */
  template<typename T>
  class ChebyshevToPowerSeries2dClosed : public ChebyshevToPowerSeries2dBase<T> {

   public:

    /// Contructor of the class.
    ChebyshevToPowerSeries2dClosed(const std::shared_ptr<ChebyshevSeries2dClosed<T>> &chebyshev) :
        ChebyshevToPowerSeries2dBase<T>(chebyshev) {
      m_x_max = chebyshev->x_max();
      m_y_max = chebyshev->y_max();
    }

   private:

    /// This method applied an affine transformation from [-1, 1] to the domain of the approximation for x.
    double AffineTransformationUnitToSegment_x(const double &xunit) const override {
      return 0.5 * (m_x_max - this->m_x_min) * xunit + 0.5 * (m_x_max + this->m_x_min);
    }

    /// This method applied an affine transformation from [-1, 1] to the domain of the approximation for y.
    double AffineTransformationUnitToSegment_y(const double &yunit) const override {
      return 0.5 * (m_y_max - this->m_y_min) * yunit + 0.5 * (m_y_max + this->m_y_min);
    }

    /// This method applied an affine transformation from the domain of the approximation to [-1, 1] for x.
    double AffineTransformationSegmentToUnit_x(const double &xdomain) const override {
      return (2. / (m_x_max - this->m_x_min)) * (xdomain - 0.5 * (m_x_max + this->m_x_min));
    }

    /// This method applied an affine transformation from the domain of the approximation to [-1, 1] for y.
    double AffineTransformationSegmentToUnit_y(const double &ydomain) const override {
      return (2. / (m_y_max - this->m_y_min)) * (ydomain - 0.5 * (m_y_max + this->m_y_min));
    }

    /// This method gives the coefficient in front of the derivate of the  double Chebychev series wrt x.
    double CoefficientDerivative_x(const double &xdomain) const {
      return 2. / (m_x_max - this->m_x_min);
    }

    /// This method gives the coefficient in front of the derivate of the  double Chebychev series wrt y.
    double CoefficientDerivative_y(const double &ydomain) const {
      return 2. / (m_y_max - this->m_y_min);
    }

   private:

    /// x maximum.
    double m_x_max;

    /// y maximum.
    double m_y_max;

  };

  /**
  * Class for computing the double power series approximation of a function from its Chebyshev series approximation over
   * [xmin, +infinity] x [ymin, +infinity].
  */
  template<typename T>
  class ChebyshevToPowerSeries2dOpened : public ChebyshevToPowerSeries2dBase<T> {

   public:

    /// Contructor of the class.
    ChebyshevToPowerSeries2dOpened(const std::shared_ptr<ChebyshevSeries2dOpened < T>> &chebyshev) :
    ChebyshevToPowerSeries2dBase <T>(chebyshev) {}

   private:

    /// This method applied an affine transformation from [-1, 1] to the domain of the approximation for x.
    double AffineTransformationUnitToSegment_x(const double &xunit) const override {
      return 2. * (this->m_x_min / (1 - xunit));
    }

    /// This method applied an affine transformation from [-1, 1] to the domain of the approximation for y.
    double AffineTransformationUnitToSegment_y(const double &yunit) const override {
      return 2. * (this->m_y_min / (1 - yunit));
    }

    /// This method applied an affine transformation from the domain of the approximation to [-1, 1] for x.
    double AffineTransformationSegmentToUnit_x(const double &xdomain) const override {
      return 1. - 2. * (this->m_x_min / xdomain);
    }

    /// This method applied an affine transformation from the domain of the approximation to [-1, 1] for y.
    double AffineTransformationSegmentToUnit_y(const double &ydomain) const override {
      return 1. - 2. * (this->m_y_min / ydomain);
    }

    /// This method gives the coefficient in front of the derivate of the  double Chebychev series wrt x.
    double CoefficientDerivative_x(const double &xdomain) const {
      return 2. * this->m_x_min / (xdomain * xdomain);
    }

    /// This method gives the coefficient in front of the derivate of the  double Chebychev series wrt y.
    double CoefficientDerivative_y(const double &ydomain) const {
      return 2. * this->m_y_min / (ydomain * ydomain);
    }

  };

}

#endif //MATHUTILS_CHEBYSHEVTOPOWERSERIES2D_H