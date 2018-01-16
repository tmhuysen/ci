//
// Created by Wulfix on 06/12/2017.
//

#include "Extras.h"

void symmatu(Eigen::MatrixXd &mat) {
    for(int x =0; x<mat.innerSize();x++){
        for (int y = x+1; y<mat.innerSize();y++){
            mat(y,x) = mat(x,y);
        }
    }


}

bool compareMat(Eigen::MatrixXd &mat1, Eigen::MatrixXd &mat2){
    if (mat1.size() != mat2.size()){
        return false;
    }else{
        for(int i = 0; i<mat1.innerSize(); i++){
            for(int j = 0; j<mat1.innerSize(); j++){
                if(abs(mat1(i,j)-mat2(i,j))>0.0001){
                    return false;
                }
            }

        }
        return true;
    }

}
