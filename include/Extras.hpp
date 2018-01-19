//
// Created by Wulfix on 06/12/2017.
//

#ifndef DOCIPROJECT_EXTRAS_H
#define DOCIPROJECT_EXTRAS_H

#include <Eigen/Dense>

struct State {
    double eigenValue;          // The energy of the solution, a.k.a. the eigenvalue
    Eigen::VectorXd eigenVector; // The coefficients of the solution with respect to the given basis, a.k.a.
    // the eigenvector corresponding to the eigenvalue
};
void symmatu(Eigen::MatrixXd &mat);
bool compareMat(Eigen::MatrixXd &mat1,Eigen::MatrixXd &mat2);
bool compareState(const State &o1, const State &o2);
bool areSame(const State &o1, const State &o2);
#endif //DOCIPROJECT_EXTRAS_H
