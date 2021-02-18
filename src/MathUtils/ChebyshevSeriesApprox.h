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
* Class for handling 2d functions for the double Chebyshev series approximation.
*/
  template<typename T>
  class Function2d {

   public:

    /// This function evaluates the function at the point (x, y).
    virtual T Evaluate(const double &x, const double &y) const = 0;

  };

/**
* Class for handling 3d functions for the double Chebyshev series approximation.
*/
  template<typename T>
  class Function3d {

   public:

    /// This function evaluates the function at the point (x, y, z).
    virtual T Evaluate(const double &x, const double &y, const double &z) const = 0;

  };

  /**
 * Class for handling and computing the double Chebyshev series approximation of a function.
 */
  template<typename T>
  class DoubleChebyshevSeriesApproxBase {

   public:

    /// Contructor of the class.
    DoubleChebyshevSeriesApproxBase(Function2d<T> *F, const double &xmin, const double &ymin, const int &order_x,
                                    const int &order_y) : m_function(F), m_x_min(xmin), m_y_min(ymin),
                                    m_order_x(order_x), m_order_y(order_y) {
      m_aij = MatrixMN<T>::Zero(order_x + 1, order_y + 1);
    }

   public:

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
        x.push_back(AffineTransformationUnitToSegment_x(tmp));
      }

      // y abscissa.
      std::vector<double> y_tilde; // In [-1,1].
      std::vector<double> y; // In [ymin,ymax].

      // y.
      for (int s = 0; s <= m_order_y; ++s) {
        double tmp = cos(MU_PI_2 * (2. * s + 1.) / (m_order_y + 1.));
        y_tilde.push_back(tmp);
        y.push_back(AffineTransformationUnitToSegment_y(tmp));
      }

      // aij.
      for (int i = 0; i <= m_order_x; ++i) {
        for (int j = 0; j <= m_order_y; ++j) {
          for (int r = 0; r <= m_order_x; ++r) {
            double xtilde = x_tilde.at(r);
            double Ti = Chebyshev_polynomial(i, xtilde);
            for (int s = 0; s <= m_order_y; ++s) {
              m_aij(i, j) += this->m_function->Evaluate(x.at(r), y.at(s)) * Ti
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

      double normal_x = AffineTransformationSegmentToUnit_x(x);
      double normal_y = AffineTransformationSegmentToUnit_y(y);

      //TODO: Vectorialiser ?
      //TODO: Etudier l'algo de Clenshaw pour des problemes de stabitilites aux bornes.
      // https://en.wikipedia.org/wiki/Clenshaw_algorithm
      // https://scicomp.stackexchange.com/questions/27865/clenshaw-type-recurrence-for-derivative-of-chebyshev-series

      T result = 0.;
      double Tn = 0.;
      double Tn_1 = 0.;
      for (int i = 0; i <= m_order_x; ++i) {

        // Definition of Ti.
        double Ti;
        if(i == 0){ // Order 0.
          Ti = 1.;
          Tn_1 = Ti;
        }
        else if(i == 1){ // Order 1.
          Ti = normal_x;
          Tn = Ti;
        }
        else{ // Order 2 and higher.
          Ti = Chebyshev_polynomial_next<double>(normal_x, Tn, Tn_1);
          Tn_1 = Tn;
          Tn = Ti;
        }

        double Tp = 0.;
        double Tp_1 = 0.;
        for (int j = 0; j <= m_order_y; ++j) {

          // Definition of Tj.
          double Tj;
          if(j == 0){ // Order 0.
            Tj = 1.;
            Tp_1 = Tj;
          }
          else if(j == 1){ // Order 1.
            Tj = normal_y;
            Tp = Tj;
          }
          else{ // Order 2 and higher.
            Tj = Chebyshev_polynomial_next<double>(normal_y, Tp, Tp_1);
            Tp_1 = Tp;
            Tp = Tj;
          }
          result += m_aij(i, j) * Ti * Tj;
        }
      }

      return result;

    }

    /// This method computes the x-derivative double Chebyshev series approximation.
    T Evaluate_derivative_x(const double &x, const double &y) const {

      double xunit = AffineTransformationSegmentToUnit_x(x);
      double yunit = AffineTransformationSegmentToUnit_y(y);

      T result = 0.;
      for (int i = 0; i <= m_order_x; ++i) {
        double Ti = Chebyshev_polynomial_derivative(i, xunit);
        for (int j = 0; j <= m_order_y; ++j) {
          result += m_aij(i, j) * Ti * Chebyshev_polynomial(j, yunit);
        }
      }
      result *= CoefficientDerivative_x(x); // With a closed segment, this coefficient is independent of normal_x.

      return result;

    }

    /// This method computes the y-derivative double Chebyshev series approximation.
    T Evaluate_derivative_y(const double &x, const double &y) const {

      double xunit = AffineTransformationSegmentToUnit_x(x);
      double yunit = AffineTransformationSegmentToUnit_y(y);

      T result = 0.;
      for (int i = 0; i <= m_order_x; ++i) {
        double Ti = Chebyshev_polynomial(i, xunit);
        for (int j = 0; j <= m_order_y; ++j) {
          result += m_aij(i, j) * Ti * Chebyshev_polynomial_derivative(j, yunit);
        }
      }
      result *= CoefficientDerivative_y(y); // With a closed segment, this coefficient is independent of normal_y.

      return result;

    }

   protected:

    /// This method applied an affine transformation from [-1, 1] to the domain of the approximation for x.
    virtual double AffineTransformationUnitToSegment_x(const double& xunit) const = 0;

    /// This method applied an affine transformation from [-1, 1] to the domain of the approximation for y.
    virtual double AffineTransformationUnitToSegment_y(const double& yunit) const = 0;

    /// This method applied an affine transformation from the domain of the approximation to [-1, 1] for x.
    virtual double AffineTransformationSegmentToUnit_x(const double& xdomain) const = 0;

    /// This method applied an affine transformation from the domain of the approximation to [-1, 1] for y.
    virtual double AffineTransformationSegmentToUnit_y(const double& ydomain) const = 0;

    /// This method gives the coefficient in front of the derivate of the  double Chebychev series wrt x.
    virtual double CoefficientDerivative_x(const double& xdomain) const = 0;

    /// This method gives the coefficient in front of the derivate of the  double Chebychev series wrt y.
    virtual double CoefficientDerivative_y(const double& ydomain) const = 0;

   protected:

    /// 2d function to be approximated.
    Function2d<T> *m_function;

    /// Order of the series approximation for x.
    int m_order_x;

    /// Order of the series approximation for y.
    int m_order_y;

    /// x minimum.
    double m_x_min;

    /// y minimum.
    double m_y_min;

    /// aij coefficients.
    MatrixMN <T> m_aij;

  };

  /**
 * Class for computing the double Chebyshev series approximation of a function over closed segments.
 */
  template<typename T>
  class DoubleChebyshevSeriesApproxClosedSegments : public DoubleChebyshevSeriesApproxBase<T> {

   public:

    /// Contructor of the class.
    DoubleChebyshevSeriesApproxClosedSegments(Function2d<T> *F, const double &xmin, const double &xmax, const double &ymin,
                          const double &ymax, const int &order_x, const int &order_y) :
                          DoubleChebyshevSeriesApproxBase<T>(F, xmin, ymin, order_x, order_y), m_x_max(xmax),
                          m_y_max(ymax) {
    }

   private:

    /// This method applied an affine transformation from [-1, 1] to the domain of the approximation for x.
    double AffineTransformationUnitToSegment_x(const double& xunit) const override {
      return 0.5 * (m_x_max - this->m_x_min) * xunit + 0.5 * (m_x_max + this->m_x_min);
    }

    /// This method applied an affine transformation from [-1, 1] to the domain of the approximation for y.
    double AffineTransformationUnitToSegment_y(const double& yunit) const override {
      return 0.5 * (m_y_max - this->m_y_min) * yunit + 0.5 * (m_y_max + this->m_y_min);
    }

    /// This method applied an affine transformation from the domain of the approximation to [-1, 1] for x.
    double AffineTransformationSegmentToUnit_x(const double& xdomain) const override {
      return (2. / (m_x_max - this->m_x_min)) * (xdomain - 0.5 * (m_x_max + this->m_x_min));
    }

    /// This method applied an affine transformation from the domain of the approximation to [-1, 1] for y.
    double AffineTransformationSegmentToUnit_y(const double& ydomain) const override {
      return (2. / (m_y_max - this->m_y_min)) * (ydomain - 0.5 * (m_y_max + this->m_y_min));
    }

    /// This method gives the coefficient in front of the derivate of the  double Chebychev series wrt x.
    double CoefficientDerivative_x(const double& xdomain) const {
      return 2. / (m_x_max - this->m_x_min);
    }

    /// This method gives the coefficient in front of the derivate of the  double Chebychev series wrt y.
    double CoefficientDerivative_y(const double& ydomain) const {
      return 2. / (m_y_max - this->m_y_min);
    }

   private:

    /// x maximum.
    double m_x_max;

    /// y maximum.
    double m_y_max;

  };

  /**
 * Class for computing the double Chebyshev series approximation of a function over half open segments.
 */
  template<typename T>
  class DoubleChebyshevSeriesApproxHalfOpenSegments : public DoubleChebyshevSeriesApproxBase<T> {

   public:

    /// Contructor of the class.
    DoubleChebyshevSeriesApproxHalfOpenSegments(Function2d<T> *F, const double &xmin, const double &ymin,
                                             const int &order_x, const int &order_y) :
        DoubleChebyshevSeriesApproxBase<T>(F, xmin, ymin, order_x, order_y) {
    }

   private:

    /// This method applied an affine transformation from [-1, 1] to the domain of the approximation for x.
    double AffineTransformationUnitToSegment_x(const double& xunit) const override {
      return 2. * (this->m_x_min / (1 - xunit));
    }

    /// This method applied an affine transformation from [-1, 1] to the domain of the approximation for y.
    double AffineTransformationUnitToSegment_y(const double& yunit) const override {
      return 2. * (this->m_y_min / (1 - yunit));
    }

    /// This method applied an affine transformation from the domain of the approximation to [-1, 1] for x.
    double AffineTransformationSegmentToUnit_x(const double& xdomain) const override {
      return 1. - 2. * (this->m_x_min / xdomain);
    }

    /// This method applied an affine transformation from the domain of the approximation to [-1, 1] for y.
    double AffineTransformationSegmentToUnit_y(const double& ydomain) const override {
      return 1. - 2. * (this->m_y_min / ydomain);
    }

    /// This method gives the coefficient in front of the derivate of the  double Chebychev series wrt x.
    double CoefficientDerivative_x(const double& xdomain) const {
      return 2. * this->m_x_min / (xdomain * xdomain);
    }

    /// This method gives the coefficient in front of the derivate of the  double Chebychev series wrt y.
    double CoefficientDerivative_y(const double& ydomain) const {
      return 2. * this->m_y_min / (ydomain * ydomain);
    }

  };

  /**
 * Class for computing the double Chebyshev series approximation of a function over a closed segment and a half open segment.
 */
  template<typename T>
  class DoubleChebyshevSeriesApproxMixedHalfOpenClosedSegments : public DoubleChebyshevSeriesApproxBase<T> {

   public:

    /// Contructor of the class.
    DoubleChebyshevSeriesApproxMixedHalfOpenClosedSegments(Function2d<T> *F, const double &xmin, const double &xmax,
                                                           const double &ymin, const int &order_x, const int &order_y) :
        DoubleChebyshevSeriesApproxBase<T>(F, xmin, ymin, order_x, order_y), m_x_max(xmax) {
    }

   private:

    /// This method applied an affine transformation from [-1, 1] to the domain of the approximation for x.
    double AffineTransformationUnitToSegment_x(const double& xunit) const override {
      return 0.5 * (m_x_max - this->m_x_min) * xunit + 0.5 * (m_x_max + this->m_x_min);
    }

    /// This method applied an affine transformation from [-1, 1] to the domain of the approximation for y.
    double AffineTransformationUnitToSegment_y(const double& yunit) const override {
      return 2. * (this->m_y_min / (1 - yunit));
    }

    /// This method applied an affine transformation from the domain of the approximation to [-1, 1] for x.
    double AffineTransformationSegmentToUnit_x(const double& xdomain) const override {
      return (2. / (m_x_max - this->m_x_min)) * (xdomain - 0.5 * (m_x_max + this->m_x_min));
    }

    /// This method applied an affine transformation from the domain of the approximation to [-1, 1] for y.
    double AffineTransformationSegmentToUnit_y(const double& ydomain) const override {
      return 1. - 2. * (this->m_y_min / ydomain);
    }

    /// This method gives the coefficient in front of the derivate of the  double Chebychev series wrt x.
    double CoefficientDerivative_x(const double& xdomain) const {
      return 2. / (m_x_max - this->m_x_min);
    }

    /// This method gives the coefficient in front of the derivate of the  double Chebychev series wrt y.
    double CoefficientDerivative_y(const double& ydomain) const {
      return 2. * this->m_y_min / (ydomain * ydomain);
    }

    /// x maximum.
    double m_x_max;

  };

  /**
* Class for handling and computing the triple Chebyshev series approximation of a function.
*/
  template<typename T>
  class TripleChebyshevSeriesApproxBase {

   public:

    /// Contructor of the class.
    TripleChebyshevSeriesApproxBase(Function3d<T> *F, const double &xmin, const double &ymin, const double &zmin,
                                    const int &order_x, const int &order_y, const int &order_z) : m_function(F),
                                    m_x_min(xmin), m_y_min(ymin), m_z_min(zmin), m_order_x(order_x), m_order_y(order_y),
                                    m_order_z(order_z) {
      m_aijk.reserve(order_z + 1);
      for (int i = 0; i <= order_z; ++i) {
        m_aijk.push_back(MatrixMN<T>::Zero(order_x + 1, order_y + 1));
      }
    }

   public:

    /// This method computes the coefficients aijk.
    void Computation_aijk() {

      // Number of coefficients in x,  y and z.
      int Nx = m_order_x + 1;
      int Ny = m_order_y + 1;
      int Nz = m_order_z + 1;

      // x abscissa.
      std::vector<double> x_tilde; // In [-1,1].
      std::vector<double> x; // In [xmin,xmax].

      // x.
      for (int r = 0; r <= m_order_x; ++r) {
        double tmp = cos(MU_PI_2 * (2. * r + 1.) / (m_order_x + 1.));
        x_tilde.push_back(tmp);
        x.push_back(AffineTransformationUnitToSegment_x(tmp));
      }

      // y abscissa.
      std::vector<double> y_tilde; // In [-1,1].
      std::vector<double> y; // In [ymin,ymax].

      // y.
      for (int s = 0; s <= m_order_y; ++s) {
        double tmp = cos(MU_PI_2 * (2. * s + 1.) / (m_order_y + 1.));
        y_tilde.push_back(tmp);
        y.push_back(AffineTransformationUnitToSegment_y(tmp));
      }

      // z abscissa.
      std::vector<double> z_tilde; // In [-1,1].
      std::vector<double> z; // In [ymin,ymax].

      // z.
      for (int t = 0; t <= m_order_z; ++t) {
        double tmp = cos(MU_PI_2 * (2. * t + 1.) / (m_order_z + 1.));
        z_tilde.push_back(tmp);
        z.push_back(AffineTransformationUnitToSegment_z(tmp));
      }

      // aijk.
      for (int i = 0; i <= m_order_x; ++i) {
        for (int j = 0; j <= m_order_y; ++j) {
          for (int k = 0; k <= m_order_z; ++k) {
            for (int r = 0; r <= m_order_x; ++r) {
              double xtilde = x_tilde.at(r);
              double Ti = Chebyshev_polynomial(i, xtilde);
              for (int s = 0; s <= m_order_y; ++s) {
                double ytilde = y_tilde.at(s);
                double Tj = Chebyshev_polynomial(j, ytilde);
                for (int t = 0; t <= m_order_z; ++t) {
                  m_aijk.at(k)(i, j) += this->m_function->Evaluate(x.at(r), y.at(s), z.at(t)) * Ti * Tj
                                        * Chebyshev_polynomial(k, z_tilde.at(t));
                }
                //TODO: Use Chebyshev_polynomial_next instead of Chebyshev_polynomial.
              }
            }
            m_aijk.at(k)(i, j) *= 8. / ((m_order_x + 1.) * (m_order_y + 1.) * (m_order_z + 1.));
            if (i == 0 and j == 0 and k == 0) {
              m_aijk.at(k)(i, j) /= 8.;
            }
            else if ((i == 0 and j != 0 and k == 0) or (i != 0 and j == 0 and k == 0)
                     or (i == 0 and j == 0 and k != 0)) {
              m_aijk.at(k)(i, j) /= 4.;
            }
            else if((i == 0 and j != 0 and k != 0) or (i != 0 and j == 0 and k != 0)
                    or (i != 0 and j != 0 and k == 0)){
              m_aijk.at(k)(i, j) /= 2.;
            }
          }
        }
      }
    }

    /// This method computes the triple Chebyshev series approximation.
    T Evaluate(const double &x, const double &y, const double &z) const {

      double normal_x = AffineTransformationSegmentToUnit_x(x);
      double normal_y = AffineTransformationSegmentToUnit_y(y);
      double normal_z = AffineTransformationSegmentToUnit_z(z);

      //TODO: Vectorialiser ?
      //TODO: Etudier l'algo de Clenshaw pour des problemes de stabitilites aux bornes.
      // https://en.wikipedia.org/wiki/Clenshaw_algorithm
      // https://scicomp.stackexchange.com/questions/27865/clenshaw-type-recurrence-for-derivative-of-chebyshev-series

      T result = 0.;
      for (int i = 0; i <= m_order_x; ++i) {
        double Ti = Chebyshev_polynomial(i, normal_x);
        for (int j = 0; j <= m_order_y; ++j) {
          double Tj = Chebyshev_polynomial(j, normal_y);
          for (int k = 0; k <= m_order_z; ++k) {
            result += m_aijk.at(k)(i, j) * Ti * Tj * Chebyshev_polynomial(k, normal_z);
          }
        }
      }

      return result;

    }

    /// This method computes the x-derivative triple Chebyshev series approximation.
    T Evaluate_derivative_x(const double &x, const double &y, const double &z) const {

      double normal_x = AffineTransformationSegmentToUnit_x(x);
      double normal_y = AffineTransformationSegmentToUnit_y(y);
      double normal_z = AffineTransformationSegmentToUnit_z(z);

      T result = 0.;
      for (int i = 0; i <= m_order_x; ++i) {
        double Ti = Chebyshev_polynomial_derivative(i, normal_x);
        for (int j = 0; j <= m_order_y; ++j) {
          double Tj = Chebyshev_polynomial(j, normal_y);
          for (int k = 0; k <= m_order_z; ++k) {
            result += m_aijk.at(k)(i, j) * Ti * Tj * Chebyshev_polynomial(k, normal_z);
          }
        }
      }
      result *= CoefficientDerivative_x(x); // With a closed segment, this coefficient is independent of normal_x.

      return result;

    }

    /// This method computes the y-derivative triple Chebyshev series approximation.
    T Evaluate_derivative_y(const double &x, const double &y, const double &z) const {

      double normal_x = AffineTransformationSegmentToUnit_x(x);
      double normal_y = AffineTransformationSegmentToUnit_y(y);
      double normal_z = AffineTransformationSegmentToUnit_z(z);

      T result = 0.;
      for (int i = 0; i <= m_order_x; ++i) {
        double Ti = Chebyshev_polynomial(i, normal_x);
        for (int j = 0; j <= m_order_y; ++j) {
          double Tj = Chebyshev_polynomial_derivative(j, normal_y);
          for (int k = 0; k <= m_order_z; ++k) {
            result += m_aijk.at(k)(i, j) * Ti * Tj * Chebyshev_polynomial(k, normal_z);
          }
        }
      }
      result *= CoefficientDerivative_y(y); // With a closed segment, this coefficient is independent of normal_x.

      return result;

    }

    /// This method computes the z-derivative triple Chebyshev series approximation.
    T Evaluate_derivative_z(const double &x, const double &y, const double &z) const {

      double normal_x = AffineTransformationSegmentToUnit_x(x);
      double normal_y = AffineTransformationSegmentToUnit_y(y);
      double normal_z = AffineTransformationSegmentToUnit_z(z);

      T result = 0.;
      for (int i = 0; i <= m_order_x; ++i) {
        double Ti = Chebyshev_polynomial(i, normal_x);
        for (int j = 0; j <= m_order_y; ++j) {
          double Tj = Chebyshev_polynomial(j, normal_y);
          for (int k = 0; k <= m_order_z; ++k) {
            result += m_aijk.at(k)(i, j) * Ti * Tj * Chebyshev_polynomial_derivative(k, normal_z);
          }
        }
      }
      result *= CoefficientDerivative_z(z); // With a closed segment, this coefficient is independent of normal_x.

      return result;

    }

   protected:

    /// This method applied an affine transformation from [-1, 1] to the domain of the approximation for x.
    virtual double AffineTransformationUnitToSegment_x(const double& xunit) const = 0;

    /// This method applied an affine transformation from [-1, 1] to the domain of the approximation for y.
    virtual double AffineTransformationUnitToSegment_y(const double& yunit) const = 0;

    /// This method applied an affine transformation from [-1, 1] to the domain of the approximation for z.
    virtual double AffineTransformationUnitToSegment_z(const double& zunit) const = 0;

    /// This method applied an affine transformation from the domain of the approximation to [-1, 1] for x.
    virtual double AffineTransformationSegmentToUnit_x(const double& xdomain) const = 0;

    /// This method applied an affine transformation from the domain of the approximation to [-1, 1] for y.
    virtual double AffineTransformationSegmentToUnit_y(const double& ydomain) const = 0;

    /// This method applied an affine transformation from the domain of the approximation to [-1, 1] for z.
    virtual double AffineTransformationSegmentToUnit_z(const double& zdomain) const = 0;

    /// This method gives the coefficient in front of the derivate of the  double Chebychev series wrt x.
    virtual double CoefficientDerivative_x(const double& xdomain) const = 0;

    /// This method gives the coefficient in front of the derivate of the  double Chebychev series wrt y.
    virtual double CoefficientDerivative_y(const double& ydomain) const = 0;

    /// This method gives the coefficient in front of the derivate of the  double Chebychev series wrt z.
    virtual double CoefficientDerivative_z(const double& zdomain) const = 0;

   protected:

    /// 3d function to be approximated.
    Function3d<T> *m_function;

    /// Order of the series approximation for x.
    int m_order_x;

    /// Order of the series approximation for y.
    int m_order_y;

    /// Order of the series approximation for z.
    int m_order_z;

    /// x minimum.
    double m_x_min;

    /// y minimum.
    double m_y_min;

    /// z minimum.
    double m_z_min;

    /// aij coefficients.
    std::vector<MatrixMN<T>> m_aijk;

  };

    /**
  * Class for computing the triple Chebyshev series approximation of a function over closed segments.
  */
  template<typename T>
  class TripleChebyshevSeriesApproxClosedSegments : public TripleChebyshevSeriesApproxBase<T> {

   public:

    /// Contructor of the class.
    TripleChebyshevSeriesApproxClosedSegments(Function3d<T> *F, const double &xmin, const double &xmax, const double &ymin,
                                const double &ymax, const double &zmin, const double &zmax, const int &order_x,
                                const int &order_y, const int &order_z) :
                                TripleChebyshevSeriesApproxBase<T>(F, xmin, ymin, zmin, order_x, order_y, order_z),
                                m_x_max(xmax), m_y_max(ymax), m_z_max(zmax) {
    }

   private:

    /// This method applied an affine transformation from [-1, 1] to the domain of the approximation for x.
    double AffineTransformationUnitToSegment_x(const double& xunit) const override {
      return 0.5 * (m_x_max - this->m_x_min) * xunit + 0.5 * (m_x_max + this->m_x_min);
    }

    /// This method applied an affine transformation from [-1, 1] to the domain of the approximation for y.
    double AffineTransformationUnitToSegment_y(const double& yunit) const override {
      return 0.5 * (m_y_max - this->m_y_min) * yunit + 0.5 * (m_y_max + this->m_y_min);
    }

    /// This method applied an affine transformation from [-1, 1] to the domain of the approximation for z.
    double AffineTransformationUnitToSegment_z(const double& zunit) const override {
      return 0.5 * (m_z_max - this->m_z_min) * zunit + 0.5 * (m_z_max + this->m_z_min);

    }

    /// This method applied an affine transformation from the domain of the approximation to [-1, 1] for x.
    double AffineTransformationSegmentToUnit_x(const double& xdomain) const override {
      return (2. / (m_x_max - this->m_x_min)) * (xdomain - 0.5 * (m_x_max + this->m_x_min));
    }

    /// This method applied an affine transformation from the domain of the approximation to [-1, 1] for y.
    double AffineTransformationSegmentToUnit_y(const double& ydomain) const override {
      return (2. / (m_y_max - this->m_y_min)) * (ydomain - 0.5 * (m_y_max + this->m_y_min));
    }

    /// This method applied an affine transformation from the domain of the approximation to [-1, 1] for z.
    double AffineTransformationSegmentToUnit_z(const double& zdomain) const override {
      return (2. / (m_z_max - this->m_z_min)) * (zdomain - 0.5 * (m_z_max + this->m_z_min));
    }

    /// This method gives the coefficient in front of the derivate of the  double Chebychev series wrt x.
    double CoefficientDerivative_x(const double& xdomain) const override {
      return 2. / (m_x_max - this->m_x_min);
    }

    /// This method gives the coefficient in front of the derivate of the  double Chebychev series wrt y.
    double CoefficientDerivative_y(const double& ydomain) const override {
      return 2. / (m_y_max - this->m_y_min);
    }

    /// This method gives the coefficient in front of the derivate of the  double Chebychev series wrt z.
    double CoefficientDerivative_z(const double& zdomain) const override {
      return 2. / (m_z_max - this->m_z_min);
    }

   private:

    /// x maximum.
    double m_x_max;

    /// y maximum.
    double m_y_max;

    /// z maximum.
    double m_z_max;

  };

  /**
* Class for computing the triple Chebyshev series approximation of a function over a closed segment in x and half-open segments in y and z.
*/
  template<typename T>
  class TripleChebyshevSeriesApproxYZHalfOpenXClosedSegments : public TripleChebyshevSeriesApproxBase<T> {

   public:

    /// Contructor of the class.
    TripleChebyshevSeriesApproxYZHalfOpenXClosedSegments(Function3d<T> *F, const double &xmin, const double &xmax,
                                                           const double &ymin, const double &zmin, const int &order_x,
                                                           const int &order_y, const int &order_z) :
        TripleChebyshevSeriesApproxBase<T>(F, xmin, ymin, zmin, order_x, order_y, order_z), m_x_max(xmax) {
    }

   private:

    /// This method applied an affine transformation from [-1, 1] to the domain of the approximation for x.
    double AffineTransformationUnitToSegment_x(const double& xunit) const override {
      return 0.5 * (m_x_max - this->m_x_min) * xunit + 0.5 * (m_x_max + this->m_x_min);
    }

    /// This method applied an affine transformation from [-1, 1] to the domain of the approximation for y.
    double AffineTransformationUnitToSegment_y(const double& yunit) const override {
      return 2. * (this->m_y_min / (1 - yunit));
    }

    /// This method applied an affine transformation from [-1, 1] to the domain of the approximation for z.
    double AffineTransformationUnitToSegment_z(const double& zunit) const override {
      return 2. * (this->m_z_min / (1 - zunit));
    }

    /// This method applied an affine transformation from the domain of the approximation to [-1, 1] for x.
    double AffineTransformationSegmentToUnit_x(const double& xdomain) const override {
      return (2. / (m_x_max - this->m_x_min)) * (xdomain - 0.5 * (m_x_max + this->m_x_min));
    }

    /// This method applied an affine transformation from the domain of the approximation to [-1, 1] for y.
    double AffineTransformationSegmentToUnit_y(const double& ydomain) const override {
      return 1. - 2. * (this->m_y_min / ydomain);
    }

    /// This method applied an affine transformation from the domain of the approximation to [-1, 1] for z.
    double AffineTransformationSegmentToUnit_z(const double& zdomain) const override {
      return 1. - 2. * (this->m_z_min / zdomain);
    }

    /// This method gives the coefficient in front of the derivate of the  double Chebychev series wrt x.
    double CoefficientDerivative_x(const double& xdomain) const override {
      return 2. / (m_x_max - this->m_x_min);
    }

    /// This method gives the coefficient in front of the derivate of the  double Chebychev series wrt y.
    double CoefficientDerivative_y(const double& ydomain) const override {
      return 2. * this->m_y_min / (ydomain * ydomain);
    }

    /// This method gives the coefficient in front of the derivate of the  double Chebychev series wrt z.
    double CoefficientDerivative_z(const double& zdomain) const override {
      return 2. * this->m_z_min / (zdomain * zdomain);
    }

   private:

    /// x maximum.
    double m_x_max;

  };

  /**
* Class for computing the triple Chebyshev series approximation of a function over a closed segment in x and half-open segments in y and z.
*/
  template<typename T>
  class TripleChebyshevSeriesApproxZHalfOpenXYClosedSegments : public TripleChebyshevSeriesApproxBase<T> {

   public:

    /// Contructor of the class.
    TripleChebyshevSeriesApproxZHalfOpenXYClosedSegments(Function3d<T> *F, const double &xmin, const double &xmax,
                                                         const double &ymin, const double &ymax, const double &zmin, const int &order_x,
                                                         const int &order_y, const int &order_z) :
        TripleChebyshevSeriesApproxBase<T>(F, xmin, ymin, zmin, order_x, order_y, order_z), m_x_max(xmax), m_y_max(ymax) {
    }

   private:

    /// This method applied an affine transformation from [-1, 1] to the domain of the approximation for x.
    double AffineTransformationUnitToSegment_x(const double& xunit) const override {
      return 0.5 * (m_x_max - this->m_x_min) * xunit + 0.5 * (m_x_max + this->m_x_min);
    }

    /// This method applied an affine transformation from [-1, 1] to the domain of the approximation for y.
    double AffineTransformationUnitToSegment_y(const double& yunit) const override {
      return 0.5 * (m_y_max - this->m_y_min) * yunit + 0.5 * (m_y_max + this->m_y_min);
    }

    /// This method applied an affine transformation from [-1, 1] to the domain of the approximation for z.
    double AffineTransformationUnitToSegment_z(const double& zunit) const override {
      return 2. * (this->m_z_min / (1 - zunit));
    }

    /// This method applied an affine transformation from the domain of the approximation to [-1, 1] for x.
    double AffineTransformationSegmentToUnit_x(const double& xdomain) const override {
      return (2. / (m_x_max - this->m_x_min)) * (xdomain - 0.5 * (m_x_max + this->m_x_min));
    }

    /// This method applied an affine transformation from the domain of the approximation to [-1, 1] for y.
    double AffineTransformationSegmentToUnit_y(const double& ydomain) const override {
      return (2. / (m_y_max - this->m_y_min)) * (ydomain - 0.5 * (m_y_max + this->m_y_min));
    }

    /// This method applied an affine transformation from the domain of the approximation to [-1, 1] for z.
    double AffineTransformationSegmentToUnit_z(const double& zdomain) const override {
      return 1. - 2. * (this->m_z_min / zdomain);
    }

    /// This method gives the coefficient in front of the derivate of the  double Chebychev series wrt x.
    double CoefficientDerivative_x(const double& xdomain) const override {
      return 2. / (m_x_max - this->m_x_min);
    }

    /// This method gives the coefficient in front of the derivate of the  double Chebychev series wrt y.
    double CoefficientDerivative_y(const double& ydomain) const override {
      return 2. / (m_y_max - this->m_y_min);
    }

    /// This method gives the coefficient in front of the derivate of the  double Chebychev series wrt z.
    double CoefficientDerivative_z(const double& zdomain) const override {
      return 2. * this->m_z_min / (zdomain * zdomain);
    }

   private:

    /// x maximum.
    double m_x_max;

    /// y maximum.
    double m_y_max;

  };

}

#endif //MATHUTILS_CHEBYSHEVSERIESAPPROX_H
