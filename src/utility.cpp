#include "utility.hpp"

void symmatu(Eigen::MatrixXd& mat) {
    for(int x =0; x<mat.innerSize();x++){
        for (int y = x+1; y<mat.innerSize();y++){
            mat(y,x) = mat(x,y);
        }
    }


}

void symmatu_reverse(Eigen::MatrixXd &mat) {
    for(int x =0; x<mat.innerSize();x++){
        for (int y = 0; y<x;y++){
            mat(y,x) = mat(x,y);
        }
    }



}
