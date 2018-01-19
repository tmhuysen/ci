//
// Created by Wulfix on 21/12/2017.
//

#define BOOST_TEST_MODULE "RHF_test"

#include "DOCI.hpp"
#include "DOCI_DENSE.hpp"
#include <hf.hpp>
#include <libwrp.hpp>
#include "DOCI_utility.hpp"

#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>



BOOST_AUTO_TEST_CASE ( DOCI_utils_test ) {
    using namespace std;
    using namespace libwrp;
    string path = "/Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/DOCILibs/DOCI_Head";
    const string xyzfilename = path + "/tests/reference_data/h2o.xyz";
    double threshold = 1.0e-06;
    string basis_name = "STO-3G";
    Molecule water (xyzfilename);
    Basis basis (water, basis_name);
    hf::rhf::RHF rhf (basis, threshold);

    CI_basis ciBasis = rhf_to_CI_basis(rhf);

}

BOOST_AUTO_TEST_CASE ( DOCI_utils2_test ) {
        using namespace std;
        using namespace libwrp;
        string path = "/Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/DOCILibs/DOCI_Head";
        const string xyzfilename = path + "/tests/reference_data/h2o.xyz";
        double threshold = 1.0e-06;
        string basis_name = "STO-3G";
        Molecule water (xyzfilename);
        Basis basis (water, basis_name);
        hf::rhf::RHF rhf (basis, threshold);

        CI_basis ciBasis = rhf_to_CI_basis(rhf);

}

BOOST_AUTO_TEST_CASE ( DOCI_DENSE_rhf_test ) {
    using namespace std;
    using namespace libwrp;
    string path = "/Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/DOCILibs/DOCI_Head";
    const string xyzfilename = path + "/tests/reference_data/h2o.xyz";
    double threshold = 1.0e-06;
    string basis_name = "STO-3G";
    Molecule water (xyzfilename);
    Basis basis (water, basis_name);
    hf::rhf::RHF rhf (basis, threshold);

    CI_basis ciBasis = rhf_to_CI_basis(rhf);

    DOCI_DENSE doci_test = DOCI_DENSE(ciBasis);
    State ground = doci_test.getGroundstates().at(0);
    double en = ground.eigenValue+ciBasis.nuc;
    cout<<endl<<en<<endl;

    BOOST_CHECK(std::abs(en - (-74.9771)) < 1.0e-04);

}