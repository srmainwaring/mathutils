//
// Created by pierre-yves on 03/03/2021.
//

#ifndef MATHUTILS_POLYNOMIAL_H
#define MATHUTILS_POLYNOMIAL_H

namespace mathutils {

  template<typename T>
  T Horner(const std::vector<double>& coefficient, const T& x0) {

    /// This method evaluates a polynomial at x0 using the Horner's method.

    int size = coefficient.size();
    T result = coefficient.at(size - 1);
    int i = size - 2;
    while (i >= 0) {
      result = result * x0 + coefficient.at(i);
      --i;
    }

    return result;

  }

  template<typename T>
  T Horner_derivative(std::vector<double>& coefficient, const T& x0) {

    /// This method evaluates the derivative of a polynomial at x0 using the Horner's method.

    int size = coefficient.size();
    T result = (size - 1) * coefficient.at(size - 1);
    int i = size - 2;
    while (i >= 1) {
      result = result * x0 + i * coefficient.at(i);
      --i;
    }

    return result;

  }

}

#endif //MATHUTILS_POLYNOMIAL_H
