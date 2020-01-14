//
// Created by pierre-yves on 14/01/20.
//

#ifndef MATHUTILS_INTEGRATION2DTRIANGLE_H
#define MATHUTILS_INTEGRATION2DTRIANGLE_H

#include "Integration2d.h"

namespace mathutils {

  /**
  * Class for computing 2d integrations over a triangle.
  */
  template<typename Tout>
  class Integration2dTriangle : public Integration2d<Tout> {

   public:

    /// Contructor of the class.
    Integration2dTriangle(Tout (*F)(Vector3d<double> x), const int &order)
        : Integration2d<Tout>(F, order) {

      // Quadrature tables.

      { // index = 0 - 1 Gauss point - order 1.
        std::vector<double> alpha = {1. / 3.};
        std::vector<double> beta = {1. / 3.};
        std::vector<double> weight = {1};
        this->m_tables->Set_table(alpha, beta, weight);
      }
      { // index = 1 - 3 Gauss points - order 2.
        std::vector<double> alpha = {0.166666666666670, 0.166666666666670, 0.666666666666670};
        std::vector<double> beta = {0.166666666666670, 0.666666666666670, 0.166666666666670};
        std::vector<double> weight = {0.333333333333330, 0.333333333333330, 0.333333333333330};
        this->m_tables->Set_table(alpha, beta, weight);
      }
      { // index = 2 - 4 Gauss point - order 3.
        std::vector<double> alpha = {0.333333333333330, 0.200000000000000, 0.200000000000000, 0.600000000000000};
        std::vector<double> beta = {0.333333333333330, 0.200000000000000, 0.600000000000000, 0.200000000000000};
        std::vector<double> weight = {-0.56250000000000, 0.520833333333330, 0.520833333333330, 0.520833333333330};
        this->m_tables->Set_table(alpha, beta, weight);
      }
      { // index = 3 - 6 Gauss points - order 4.
        std::vector<double> alpha = {0.445948490915970, 0.445948490915970, 0.108103018168070, 0.091576213509770,
                                     0.091576213509770, 0.816847572980460};
        std::vector<double> beta = {0.445948490915970, 0.108103018168070, 0.445948490915970, 0.091576213509770,
                                    0.816847572980460, 0.091576213509770};
        std::vector<double> weight = {0.223381589678010, 0.223381589678010, 0.223381589678010, 0.109951743655320,
                                      0.109951743655320, 0.109951743655320};
        this->m_tables->Set_table(alpha, beta, weight);
      }
      { // index = 4 - 7 Gauss points - order 5.
        std::vector<double> alpha = {0.333333333333330, 0.470142064105110, 0.470142064105110, 0.059715871789770,
                                     0.101286507323460, 0.101286507323460, 0.797426985353090};
        std::vector<double> beta = {0.333333333333330, 0.470142064105110, 0.059715871789770, 0.470142064105110,
                                    0.101286507323460, 0.797426985353090, 0.101286507323460};
        std::vector<double> weight = {0.225000000000000, 0.132394152788510, 0.132394152788510, 0.132394152788510,
                                      0.125939180544830, 0.125939180544830, 0.125939180544830};
        this->m_tables->Set_table(alpha, beta, weight);
      }

      // Warning.
      if (this->m_order > this->m_tables->get_size_tables()) {
        std::cout << "\nToo high order of quadrature for triangle. Use lower order." << std::endl;
        exit(0);
      }
      if (this->m_order < 1) {
        std::cout << "\nToo low order of quadrature for triangle. Use higher order." << std::endl;
        exit(0);
      }

    }

    /// This function computes the 2d integration over a triangle.
    Tout Compute(const Vector3d<double> &vertex_1, const Vector3d<double> &vertex_2,
                 const Vector3d<double> &vertex_3) override {

      // Face area.
      Vector3d<double> e1 = vertex_1 - vertex_3;
      Vector3d<double> e2 = vertex_2 - vertex_3;
      double area = 0.5 * ((e1.cross(e2))).norm();

      // Computation of the 2d integration over a standart triangle.
      Tout integral = ComputeOverStandartTriangle(vertex_1, vertex_2, vertex_3);

      // 2d integration over the real triangle.
      integral *= 2*area;

      return(integral);

    }

   private:

    /// This function computes the 2d integration over a standart triangle.
    Tout ComputeOverStandartTriangle(const Vector3d<double> &vertex_1, const Vector3d<double> &vertex_2,
                                     const Vector3d<double> &vertex_3) {

      // Quadrature coefficients.
      std::vector<double> alpha = this->m_tables->alpha(this->m_order);
      std::vector<double> beta = this->m_tables->beta(this->m_order);
      std::vector<double> weight = this->m_tables->weight(this->m_order);

      // Computation of the 2d integration.
      Tout result = 0;
      Tout tmp_val;
      for (unsigned int i = 0; i < weight.size(); i++) {

        // Position of the Gauss points.
        Vector3d<double> x = (1 - alpha.at(i) - beta.at(i)) * vertex_1 + alpha.at(i) * vertex_2 + beta.at(i) * vertex_3;

        tmp_val = (*this->m_integrand)(x);
        tmp_val *= weight.at(i) * 0.5; // The 1/2 coefficient matches the area of the standart triangle.
        result += tmp_val;
      }

      return (result);

    }

  };

}

#endif //MATHUTILS_INTEGRATION2DTRIANGLE_H
