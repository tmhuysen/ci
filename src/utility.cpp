#include "utility.hpp"

void doci::symmatu(Eigen::MatrixXd& mat) {
    for(int x =0; x<mat.innerSize();x++){
        for (int y = x+1; y<mat.innerSize();y++){
            mat(y,x) = mat(x,y);
        }
    }
}

void doci::symmatu_reverse(Eigen::MatrixXd &mat) {
    for(int x =0; x<mat.innerSize();x++){
        for (int y = 0; y<x;y++){
            mat(y,x) = mat(x,y);
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
/*
// trim from both ends (in place)
static inline void doci::trim(std::string &s) {
    doci::ltrim(s);
    doci::rtrim(s);
}

// trim from start (copying)
static inline std::string doci::ltrim_copy(std::string s) {
    doci::ltrim(s);
    return s;
}

// trim from end (copying)
static inline std::string doci::rtrim_copy(std::string s) {
    doci::rtrim(s);
    return s;
}

// trim from both ends (copying)
static inline std::string trim_copy(std::string s) {
    doci::trim(s);
    return s;
}
 */