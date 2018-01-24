#define BOOST_TEST_MODULE "CI_basis_test"

#include "CI_basis.hpp"

#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>


BOOST_AUTO_TEST_CASE ( constructor ) {
    // FIXME: This test should be finished
//    using namespace std;
//    using namespace libwint;
//    string path = "/Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/DOCILibs/DOCI_Head";
//    const string doci = path + "/tests/reference_data/doci_ref/h2o.xyz";
//
//
//
//    CI_basis ciBasis = file_to_CI_basis(doci);

}


BOOST_AUTO_TEST_CASE ( constructor_rhf ) {
    const std::string xyzfilename = "../tests/reference_data/h2o.xyz";
    double threshold = 1.0e-06;
    std::string basis_name = "STO-3G";
    libwint::Molecule water (xyzfilename);
    libwint::Basis basis (water, basis_name);
    hf::rhf::RHF rhf (basis, threshold);

    doci::CI_basis ciBasis (rhf);
}


BOOST_AUTO_TEST_CASE ( constructor_filename) {
    const std::string doci = "../tests/reference_data/doci_ref/beh_cation_ref.txt";
    doci::CI_basis ciBasis(doci);

}