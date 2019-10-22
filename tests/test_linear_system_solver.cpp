//
// Created by pierre-yves on 22/10/19.
//

#include <iostream>
#include "MathUtils/MathUtils.h"

#include <iostream>
#include <Eigen/Dense>
using namespace std;
using namespace Eigen;


int main(int argc, char* argv[]) {

  // Matrix66.
  mathutils::Matrix66<double> A66;
  A66 << 2,5,7,9,1,-6,
       3,7,-4,8,-6,1,
       -10,13,24,-14,3,0,
       3,0,4,8,9,-6,
       1,7,9,3,-20,4,
       0,0,3,4,8,9;

  mathutils::Vector6d<double> b6;
  b6 << 1, 2, 3, 4, 5, 6;

  auto x6 = A66.LUSolver<mathutils::Vector6d<double>, mathutils::Vector6d<double>>(b6);

  std::cout << "" << std::endl;
  std::cout << x6 << std::endl;

  // Matrix33.
  mathutils::Matrix33<double> A33;
  A33 << 2,5,7,
      3,7,-4,
      -10,13,24;

  mathutils::Vector3d<double> b3;
  b3 << 1, 2, 3;

  auto x3 = A33.LUSolver<mathutils::Vector3d<double>, mathutils::Vector3d<double>>(b3);

  std::cout << "" << std::endl;
  std::cout << x3 << std::endl;

  // MatrixMN.
  mathutils::MatrixMN<double> AMN = mathutils::MatrixMN<double>::Zero(6, 6);
  AMN = A66;

  mathutils::MatrixMN<double> bN = mathutils::MatrixMN<double>::Zero(6, 2);
  bN << 1, 1, 2, 2, 3, 3,
      4, 4, 5, 5, 6, 6;

  auto xN = AMN.LUSolver<mathutils::MatrixMN<double>, mathutils::MatrixMN<double>>(bN);

  std::cout << "" << std::endl;
  std::cout << xN << std::endl;

  return 0;
}