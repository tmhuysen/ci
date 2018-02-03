
#include "utility.hpp"

void doci::symmatu(Eigen::MatrixXd& mat) { // no longer used
    for(int x =0; x<mat.innerSize();x++){
        for (int y = x+1; y<mat.innerSize();y++){
            mat(y,x) = mat(x,y);
        }
    }
}

void doci::symmatu_reverse(Eigen::MatrixXd &mat) { // no longer used
    for(int x =0; x<mat.innerSize();x++){
        for (int y = 0; y<x;y++){
            mat(y,x) = mat(x,y);
        }
    }
}

void doci::symmatu_tensor_reverse(Eigen::Tensor<double, 4>& tei) { //FIXME doesn't work (only for sparse DOCI tensors
    for (size_t i = 0; i<tei.dimension(0); i++) {
        for (size_t j = 0; j <= i; j++) {
            for (size_t k = 0; k <= i; k++) {
                for (size_t l = 0; l <= k && l<= j &&10*k +l<= 10*i+j; l++) {
                    tei(j, i, k, l) = tei(i, j, k, l);
                    tei(j, i, l, k) = tei(i, j, k, l);
                    tei(i, j, l, k) = tei(i, j, k, l);

                }
            }
        }
    }


}

// trim from start (in place)
void doci::ltrim(std::string &s) {
    s.erase(s.begin(), std::find_if(s.begin(), s.end(),
                                    std::not1(std::ptr_fun<int, int>(std::isspace))));
}

// trim from end (in place)
void doci::rtrim(std::string &s) {
    s.erase(std::find_if(s.rbegin(), s.rend(),
                         std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
}

