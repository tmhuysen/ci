#define BOOST_TEST_MODULE "Utility_tests"

#include "utility.hpp"

#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>




BOOST_AUTO_TEST_CASE ( symmatu_tensor_4) {
    //IGNORE
    int K = 3;
    Eigen::Tensor<double, 4> test_tensor(K,K,K,K);
    test_tensor.setZero();
    for (size_t i = 0; i<test_tensor.dimension(0); i++) {
        for (size_t j = 0; j <= i; j++) {
            for (size_t k = 0; k <= i; k++) {
                for (size_t l = 0; l <= k && l<= j &&10*k +l<= 10*i+j; l++) {
                    size_t a =i;size_t b=j;size_t c=k;size_t d=l;
                    //std::cout<<std::endl<<" "<<i<<" "<<j<<" "<<k<<" "<<l;
                    test_tensor(a,b,c,d) = 1;
                    test_tensor(b, a, c, d) = 1; //test_tensor(i, j, k, l);
                    test_tensor(b, a, d, c) = 1;//test_tensor(i, j, k, l);
                    test_tensor(a, b, d, c) = 1;//test_tensor(i, j, k, l);

                    a = l;b=k;c=j;d=i;
                    test_tensor(a,b,c,d) = 1;
                    test_tensor(b, a, c, d) = 1; //test_tensor(i, j, k, l);
                    test_tensor(b, a, d, c) = 1;//test_tensor(i, j, k, l);
                    test_tensor(a, b, d, c) = 1;//test_tensor(i, j, k, l);

                    a = l;b=k;c=i;d=j;
                    test_tensor(a,b,c,d) = 1;
                    test_tensor(b, a, c, d) = 1; //test_tensor(i, j, k, l);
                    test_tensor(b, a, d, c) = 1;//test_tensor(i, j, k, l);
                    test_tensor(a, b, d, c) = 1;//test_tensor(i, j, k, l);

                    a = l;b=i;c=j;d=k;
                    test_tensor(a,b,c,d) = 1;
                    test_tensor(b, a, c, d) = 1; //test_tensor(i, j, k, l);
                    test_tensor(b, a, d, c) = 1;//test_tensor(i, j, k, l);
                    test_tensor(a, b, d, c) = 1;//test_tensor(i, j, k, l);


                }
            }
        }
    }

    for (size_t i = 0; i<test_tensor.dimension(0); i++) {
        for (size_t j = 0; j < test_tensor.dimension(0); j++) {
            for (size_t k = 0; k <test_tensor.dimension(0); k++) {
                for (size_t l = 0; l <test_tensor.dimension(0); l++) {
                    if(!test_tensor(i,j,k,l)){
                        //std::cout<<std::endl<<" "<<i<<" "<<j<<" "<<k<<" "<<l;

                    }



                }
            }
        }
    }


}