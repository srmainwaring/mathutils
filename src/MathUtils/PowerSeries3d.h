//
// Created by pierre-yves on 04/03/2021.
//

#ifndef MATHUTILS_POWERSERIES3D_H
#define MATHUTILS_POWERSERIES3D_H

#include "Matrix.h"
#include "Polynomial.h"
#include "Functions.h"

namespace mathutils {

/**
* Class for handling and computing the triple power series approximation of a function.
*/
  template<typename T>
  class PowerSeries3dBase {

   public:

    /// Contructor of the class.
    PowerSeries3dBase(const std::vector<MatrixMN <T>> &bijk, const double &xmin, const double &ymin, const double &zmin)
    : m_bijk(bijk), m_x_min(xmin), m_y_min(ymin), m_z_min(zmin) {

      m_order_x = bijk.at(0).rows() - 1;
      m_order_y = bijk.at(0).cols() - 1;
      m_order_z = bijk.size() - 1;
    }

    /// Getter for the series coefficients.
    std::vector<MatrixMN <T>> bijk() {
      return m_bijk;
    }

    /// Getter of m_x_min.
    double xmin() {
      return m_x_min;
    }

    /// Getter of m_y_min.
    double ymin() {
      return m_y_min;
    }

    /// Getter of m_z_min.
    double zmin() {
      return m_z_min;
    }

    /// Getter of m_order_x.
    int order_x() {
      return m_order_x;
    }

    /// Getter of m_order_y.
    int order_y() {
      return m_order_y;
    }

    /// Getter of m_order_z.
    int order_z() {
      return m_order_z;
    }

    /// Getter of m_order_z.
    std::string type() {
      return m_type;
    }

    /// This method computes the triple power series approximation for unit coordinates.
    T Evaluate_unit(const double &xunit, const double &yunit, const double &zunit) const {

      // Partial sum.
      std::vector<double> rk;
      rk.reserve(m_order_z + 1);
      for (int k = 0; k <= m_order_z; ++k) {
        std::vector<double> qi;
        qi.reserve(m_order_x + 1);
        for (int i = 0; i <= m_order_x; ++i) {
          Eigen::VectorXd tmp = m_bijk.at(k).row(i);
          std::vector<double> pj(tmp.data(), tmp.data() + tmp.size());
          qi.push_back(Horner<double>(pj, yunit));
        }
        rk.push_back(Horner<double>(qi, xunit));
      }
      T result = Horner<double>(rk, zunit);

      return result;

    }

    /// This method computes the triple power series approximation.
    T Evaluate(const double &x, const double &y, const double &z) const {

      // Parameters.
      double xunit = AffineTransformationSegmentToUnit_x(x);
      double yunit = AffineTransformationSegmentToUnit_y(y);
      double zunit = AffineTransformationSegmentToUnit_z(z);

      return Evaluate_unit(xunit, yunit, zunit);

    }

    /// This method computes the x-derivative triple power series approximation for unit coordinates.
    T Evaluate_derivative_x_unit(const double &x, const double &xunit, const double &yunit, const double &zunit) const {

      // Partial sum.
      std::vector<double> rk;
      rk.reserve(m_order_z + 1);
      for (int k = 0; k <= m_order_z; ++k) {
        std::vector<double> qi;
        qi.reserve(m_order_x + 1);
        for (int i = 0; i <= m_order_x; ++i) {
          Eigen::VectorXd tmp = m_bijk.at(k).row(i);
          std::vector<double> pj(tmp.data(), tmp.data() + tmp.size());
          qi.push_back(Horner<double>(pj, yunit));
        }
        rk.push_back(Horner_derivative<double>(qi, xunit));
      }
      T result = CoefficientDerivative_x(x) * Horner<double>(rk, zunit);

      return result;

    }

    /// This method computes the x-derivative triple power series approximation.
    T Evaluate_derivative_x(const double &x, const double &y, const double &z) const {

      // Parameters.
      double xunit = AffineTransformationSegmentToUnit_x(x);
      double yunit = AffineTransformationSegmentToUnit_y(y);
      double zunit = AffineTransformationSegmentToUnit_z(z);

      return Evaluate_derivative_x_unit(x, xunit, yunit, zunit);

    }

    /// This method computes the y-derivative triple power series approximation for unit coordinates.
    T Evaluate_derivative_y_unit(const double &y, const double &xunit, const double &yunit, const double &zunit) const {

      // Partial sum.
      std::vector<double> rk;
      rk.reserve(m_order_z + 1);
      for (int k = 0; k <= m_order_z; ++k) {
        std::vector<double> qi;
        qi.reserve(m_order_x + 1);
        for (int i = 0; i <= m_order_x; ++i) {
          Eigen::VectorXd tmp = m_bijk.at(k).row(i);
          std::vector<double> pj(tmp.data(), tmp.data() + tmp.size());
          qi.push_back(Horner_derivative<double>(pj, yunit));
        }
        rk.push_back(Horner<double>(qi, xunit));
      }
      T result = CoefficientDerivative_y(y) * Horner<double>(rk, zunit);

      return result;

    }

    /// This method computes the y-derivative triple power series approximation.
    T Evaluate_derivative_y(const double &x, const double &y, const double &z) const {

      // Parameters.
      double xunit = AffineTransformationSegmentToUnit_x(x);
      double yunit = AffineTransformationSegmentToUnit_y(y);
      double zunit = AffineTransformationSegmentToUnit_z(z);

      return Evaluate_derivative_y_unit(y, xunit, yunit, zunit);

    }

    /// This method computes the z-derivative triple power series approximation for unit coordinates.
    T Evaluate_derivative_z_unit(const double &z, const double &xunit, const double &yunit, const double &zunit) const {

      // Partial sum.
      std::vector<double> rk;
      rk.reserve(m_order_z + 1);
      for (int k = 0; k <= m_order_z; ++k) {
        std::vector<double> qi;
        qi.reserve(m_order_x + 1);
        for (int i = 0; i <= m_order_x; ++i) {
          Eigen::VectorXd tmp = m_bijk.at(k).row(i);
          std::vector<double> pj(tmp.data(), tmp.data() + tmp.size());
          qi.push_back(Horner<double>(pj, yunit));
        }
        rk.push_back(Horner<double>(qi, xunit));
      }
      T result = CoefficientDerivative_z(z) * Horner_derivative<double>(rk, zunit);

      return result;

    }

    /// This method computes the z-derivative triple power series approximation.
    T Evaluate_derivative_z(const double &x, const double &y, const double &z) const {

      // Parameters.
      double xunit = AffineTransformationSegmentToUnit_x(x);
      double yunit = AffineTransformationSegmentToUnit_y(y);
      double zunit = AffineTransformationSegmentToUnit_z(z);

      return Evaluate_derivative_z_unit(z, xunit, yunit, zunit);

    }

   protected:

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

    /// Pij coefficients.
    std::vector<MatrixMN<T>> m_bijk;

    /// Type of approximation.
    std::string m_type;

  };

  /**
 * Class for computing the triple power series approximation of a function over [xmin, xmax] x [ymin, ymax] x [zmin, zmax].
 */
  template<typename T>
  class PowerSeries3dClosed : public PowerSeries3dBase<T> {

   public:

    /// Contructor of the class.
    PowerSeries3dClosed(const std::vector<MatrixMN <T>> &bijk, const double &xmin, const double &xmax, const double &ymin,
                        const double &ymax, const double &zmin, const double &zmax) :
        PowerSeries3dBase<T>(bijk, xmin, ymin, zmin), m_x_max(xmax), m_y_max(ymax), m_z_max(zmax) {
      this->m_type = "PowerSeries3dClosed";
    }

    /// Getter of m_x_max.
    double xmax() {
      return m_x_max;
    }

    /// Getter of m_y_max.
    double ymax() {
      return m_y_max;
    }

    /// Getter of m_z_max.
    double zmax() {
      return m_z_max;
    }

   private:

    /// This method applied an affine transformation from the domain of the approximation to [-1, 1] for x.
    double AffineTransformationSegmentToUnit_x(const double &xdomain) const override {
      return (2. / (m_x_max - this->m_x_min)) * (xdomain - 0.5 * (m_x_max + this->m_x_min));
    }

    /// This method applied an affine transformation from the domain of the approximation to [-1, 1] for y.
    double AffineTransformationSegmentToUnit_y(const double &ydomain) const override {
      return (2. / (m_y_max - this->m_y_min)) * (ydomain - 0.5 * (m_y_max + this->m_y_min));
    }

    /// This method applied an affine transformation from the domain of the approximation to [-1, 1] for z.
    double AffineTransformationSegmentToUnit_z(const double &zdomain) const override {
      return (2. / (m_z_max - this->m_z_min)) * (zdomain - 0.5 * (m_z_max + this->m_z_min));
    }

    /// This method gives the coefficient in front of the derivate of the  double Chebychev series wrt x.
    double CoefficientDerivative_x(const double &xdomain) const override {
      return 2. / (m_x_max - this->m_x_min);
    }

    /// This method gives the coefficient in front of the derivate of the  double Chebychev series wrt y.
    double CoefficientDerivative_y(const double &ydomain) const override {
      return 2. / (m_y_max - this->m_y_min);
    }

    /// This method gives the coefficient in front of the derivate of the  double Chebychev series wrt z.
    double CoefficientDerivative_z(const double &zdomain) const override {
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
* Class for computing the triple power series approximation of a function over [xmin, xmax] x [ymin, +infinity] x [zmin, +infinity].
*/
  template<typename T>
  class PowerSeries3dYZOpenedXClosed : public PowerSeries3dBase<T> {

   public:

    /// Contructor of the class.
    PowerSeries3dYZOpenedXClosed(const std::vector<MatrixMN <T>> &bijk, const double &xmin, const double &xmax,
                                 const double &ymin, const double &zmin) :
    PowerSeries3dBase<T>(bijk, xmin, ymin, zmin), m_x_max(xmax) {
      this->m_type = "PowerSeries3dYZOpenedXClosed";
    }

    /// Getter of m_x_max.
    double xmax() {
      return m_x_max;
    }

   private:

    /// This method applied an affine transformation from the domain of the approximation to [-1, 1] for x.
    double AffineTransformationSegmentToUnit_x(const double &xdomain) const override {
      return (2. / (m_x_max - this->m_x_min)) * (xdomain - 0.5 * (m_x_max + this->m_x_min));
    }

    /// This method applied an affine transformation from the domain of the approximation to [-1, 1] for y.
    double AffineTransformationSegmentToUnit_y(const double &ydomain) const override {
      return 1. - 2. * (this->m_y_min / ydomain);
    }

    /// This method applied an affine transformation from the domain of the approximation to [-1, 1] for z.
    double AffineTransformationSegmentToUnit_z(const double &zdomain) const override {
      return 1. - 2. * (this->m_z_min / zdomain);
    }

    /// This method gives the coefficient in front of the derivate of the  double Chebychev series wrt x.
    double CoefficientDerivative_x(const double &xdomain) const override {
      return 2. / (m_x_max - this->m_x_min);
    }

    /// This method gives the coefficient in front of the derivate of the  double Chebychev series wrt y.
    double CoefficientDerivative_y(const double &ydomain) const override {
      return 2. * this->m_y_min / (ydomain * ydomain);
    }

    /// This method gives the coefficient in front of the derivate of the  double Chebychev series wrt z.
    double CoefficientDerivative_z(const double &zdomain) const override {
      return 2. * this->m_z_min / (zdomain * zdomain);
    }

   private:

    /// x maximum.
    double m_x_max;

  };

  /**
* Class for computing the triple power series approximation of a function over [xmin, xmax] x [ymin, ymx] x [zmin, +infinity].
*/
  template<typename T>
  class PowerSeries3dZOpenedXYClosed : public PowerSeries3dBase<T> {

   public:

    /// Contructor of the class.
    PowerSeries3dZOpenedXYClosed(const std::vector<MatrixMN <T>> &bijk, const double &xmin, const double &xmax,
                                 const double &ymin, const double &ymax, const double &zmin) :
    PowerSeries3dBase<T>(bijk, xmin, ymin, zmin), m_x_max(xmax), m_y_max(ymax) {
      this->m_type = "PowerSeries3dZOpenedXYClosed";
    }

    /// Getter of m_x_max.
    double xmax() {
      return m_x_max;
    }

    /// Getter of m_y_max.
    double ymax() {
      return m_y_max;
    }

   private:

    /// This method applied an affine transformation from the domain of the approximation to [-1, 1] for x.
    double AffineTransformationSegmentToUnit_x(const double &xdomain) const override {
      return (2. / (m_x_max - this->m_x_min)) * (xdomain - 0.5 * (m_x_max + this->m_x_min));
    }

    /// This method applied an affine transformation from the domain of the approximation to [-1, 1] for y.
    double AffineTransformationSegmentToUnit_y(const double &ydomain) const override {
      return (2. / (m_y_max - this->m_y_min)) * (ydomain - 0.5 * (m_y_max + this->m_y_min));
    }

    /// This method applied an affine transformation from the domain of the approximation to [-1, 1] for z.
    double AffineTransformationSegmentToUnit_z(const double &zdomain) const override {
      return 1. - 2. * (this->m_z_min / zdomain);
    }

    /// This method gives the coefficient in front of the derivate of the  double Chebychev series wrt x.
    double CoefficientDerivative_x(const double &xdomain) const override {
      return 2. / (m_x_max - this->m_x_min);
    }

    /// This method gives the coefficient in front of the derivate of the  double Chebychev series wrt y.
    double CoefficientDerivative_y(const double &ydomain) const override {
      return 2. / (m_y_max - this->m_y_min);
    }

    /// This method gives the coefficient in front of the derivate of the  double Chebychev series wrt z.
    double CoefficientDerivative_z(const double &zdomain) const override {
      return 2. * this->m_z_min / (zdomain * zdomain);
    }

   private:

    /// x maximum.
    double m_x_max;

    /// y maximum.
    double m_y_max;

  };

}

#endif //MATHUTILS_POWERSERIES3D_H
