#include "CI_basis.hpp"


/** Default constructor
 */
doci::CI_basis::CI_basis() {}


/** Constructor based on a given RHF instance
 */
doci::CI_basis::CI_basis(hf::rhf::RHF& rhf) {
    // Transform the one- and two-electron integrals from the AO basis to the SO basis
    Eigen::MatrixXd SO_V = libwint::transform_AO_to_SO(rhf.basis.V,rhf.C_canonical);
    Eigen::MatrixXd SO_T = libwint::transform_AO_to_SO(rhf.basis.T,rhf.C_canonical);
    this->one_ints = SO_V + SO_T;
    this->two_ints = libwint::transform_AO_to_SO(rhf.basis.tei,rhf.C_canonical);

    // Set some more parameters
    this->internuclear_repulsion = rhf.basis.molecule.internuclear_repulsion();
    this->K = rhf.basis.nbf();
    this->nelec = rhf.basis.molecule.nelec;
}


/** Constructor based on a given filename
 *  Very specific file.
 */
doci::CI_basis::CI_basis(const std::string &filename) {
    std::ifstream is (filename);

    if (is.is_open()) {
        std::string startline; // first line contains orbitals and electron count.
        std::getline(is, startline); // extract the line from the file.
        std::stringstream linestream(startline); // turn the string into a stream.

        int value;
        char itter;
        int counter = 0;
        while(linestream >> itter && counter<2){
            if(itter == '='){ //
                if(counter == 0){
                    linestream >> value;
                    this->K = value;
                    counter++;
                }else{
                    linestream >> value;
                    this->nelec = value;
                    counter++;
                }
            }

        }

        this->two_ints = Eigen::Tensor<double, 4>(K,K,K,K);
        this->two_ints.setZero();
        this->one_ints = Eigen::MatrixXd::Zero(K,K);

        for(int i =0;i<3;i++){
            std::getline(is, startline); //skip 3 lines;

        }

        std::string intline;


        while(std::getline(is, intline)){
            std::stringstream intss(intline);
            std::vector<double> inputs;
            double test;
            while(intss >> test){
                inputs.push_back(test);

            }

            double integral = inputs.at(0);
            size_t index1 = inputs.at(1);
            size_t index2 = inputs.at(2);
            size_t index3 = inputs.at(3);
            size_t index4 = inputs.at(4);


            if(index1>0 && index2>0  && index3>0  && index4>0 ){
                size_t i = index1-1;
                size_t j = index2-1;
                size_t l = index3-1;
                size_t k = index4-1;
                this->two_ints(i,j,k,l) = integral;
                this->two_ints(j, i, k, l) = this->two_ints(i, j, k, l);
                this->two_ints(j, i, l, k) = this->two_ints(i, j, k, l);
                this->two_ints(i, j, l, k) = this->two_ints(i, j, k, l);
                //FIXME only works for DOCI skips many indexes.

            }
            if(index1>0 && index2>0 && index3==0){
                size_t i = index1-1;
                size_t j = index2-1;
                this->one_ints(i,j) = integral;


            }
            if(index1==0){
                this->internuclear_repulsion = integral;

            }

        }

    }else{
        std::cout<<"can't open file";
    }

}


/**
 * Getters
 */
double doci::CI_basis::getOne_int(size_t index1, size_t index2) const {
    return one_ints(index1,index2);
}

double doci::CI_basis::getTwo_int(size_t index1, size_t index2, size_t index3, size_t index4) const {
    return two_ints(index1,index2,index3,index4);
}

double doci::CI_basis::getInternuclear_repulsion() const {
    return internuclear_repulsion;
}

size_t doci::CI_basis::getK() const {
    return K;
}

size_t doci::CI_basis::getNelec() const {
    return nelec;
}
