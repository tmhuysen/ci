#define BOOST_TEST_MODULE "RHF_test"

#include "DOCI.hpp"
#include "DOCI_DENSE.hpp"
#include <hf.hpp>
#include <libwrp.hpp>
#include "DOCI_utility.hpp"

#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>



BOOST_AUTO_TEST_CASE ( file_utils_test ) {
    using namespace std;
    using namespace libwrp;
    string path = "/Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/DOCILibs/DOCI_Head";
    const string doci = path + "/tests/reference_data/doci_ref/h2o.xyz";



    CI_basis ciBasis = file_to_CI_basis(doci);

}

