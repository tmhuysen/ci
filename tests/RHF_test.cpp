//
// Created by Wulfix on 21/12/2017.
//

#define BOOST_TEST_MODULE "RHF_test"

#include "DOCI.h"
#include "RHFWrapper.h"

#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>
BOOST_AUTO_TEST_CASE ( RHF_int_test ) {
    using namespace std;
    using namespace libwrp;
    string path = "/Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/DOCILibs/DOCI_Head";
    const string xyzfilename = path + "/tests/reference_data/h2o.xyz";
    double threshold = 1.0e-06;
    string basis_name = "STO-3G";
    Molecule water (xyzfilename);
    Basis basis (water, basis_name);
    hf::rhf::RHF rhf (basis, threshold);
    double ref_energy = -74.942080055631;
    BOOST_CHECK(std::abs(rhf.energy - ref_energy) < 1.0e-06);






}