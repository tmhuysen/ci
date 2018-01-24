#define BOOST_TEST_MODULE "DenseDOCI_test"

#include "DenseDOCI.hpp"

#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>



BOOST_AUTO_TEST_CASE ( first_test ) {  // FIXME use better name
    const std::string xyzfilename = "../tests/reference_data/h2o.xyz";
    double threshold = 1.0e-06;
    std::string basis_name = "STO-3G";
    libwint::Molecule water (xyzfilename);
    libwint::Basis basis (water, basis_name);
    hf::rhf::RHF rhf (basis, threshold);

    doci::CI_basis ciBasis (rhf);

    doci::DenseDOCI doci_test (ciBasis);
    doci::State ground = doci_test.getGroundstates().at(0);
    double en = ground.eval + ciBasis.getInternuclear_repulsion();
    std::cout << std::endl << en << std::endl;

    BOOST_CHECK(std::abs(en - (-74.9771)) < 1.0e-04);  // FIXME: add reference to DOCI energy in comments
}

BOOST_AUTO_TEST_CASE ( second_test ) {  // FIXME use better name
    const std::string doci = "../tests/reference_data/doci_ref/beh_cation_ref.txt";
    doci::CI_basis ciBasis(doci);

    doci::DenseDOCI doci_test (ciBasis);
    doci::State ground = doci_test.getGroundstates().at(0);
    double en = ground.eval + ciBasis.getInternuclear_repulsion();
    std::cout << std::endl << en << std::endl;
    
    //BOOST_CHECK(std::abs(en - (-74.9771)) < 1.0e-04);  // FIXME: add reference to DOCI energy in comments
}