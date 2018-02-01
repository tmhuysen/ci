#define BOOST_TEST_MODULE "SparseDOCI_test"

#include "SparseDOCI.hpp"
#include "utility.hpp"

#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>



BOOST_AUTO_TEST_CASE ( self_test_sparse ) {
    const std::string xyzfilename = "../tests/reference_data/h2o.xyz";
    double threshold = 1.0e-06;
    std::string basis_name = "STO-3G";
    libwint::Molecule water (xyzfilename);
    libwint::Basis basis (water, basis_name);
    hf::rhf::RHF rhf (basis, threshold);
    doci::CI_basis ciBasis (rhf);

    doci::SparseDOCI doci_test (ciBasis);
    doci::State ground = doci_test.getGroundstates().at(0);
    double en = ground.getEval() + ciBasis.getInternuclear_repulsion();

    BOOST_CHECK(std::abs(en - (-74.9771)) < 1.0e-04);  // compares to energy with the program's own previous calc.
}


BOOST_AUTO_TEST_CASE ( doci_ref_beh_test_sparse ) {
    //beh_cation
    const std::string doci = "../tests/reference_data/doci_ref/beh_cation_ap1rog_631g.txt";
    //loading the basis
    doci::CI_basis ciBasis(doci);
    //calculating doci
    doci::SparseDOCI doci_test(ciBasis);
    //getting the electronic ground energy
    doci::State ground = doci_test.getGroundstates().at(0);
    //getting the total energy
    double en = ground.getEval() + ciBasis.getInternuclear_repulsion();
    //compare
    BOOST_CHECK(std::abs(en - (-14.8782216937)) < 1.0e-06); //energy from beh_cation
}

BOOST_AUTO_TEST_CASE ( doci_ref_lih_test_sparse ) {
    //lih
    const std::string doci = "../tests/reference_data/doci_ref/lih_ap1rog_631g.txt";
    //loading the basis
    doci::CI_basis ciBasis(doci);
    //calculating doci
    doci::SparseDOCI doci_test(ciBasis);
    //getting the electronic ground energy
    doci::State ground = doci_test.getGroundstates().at(0);
    //getting the total energy
    double en = ground.getEval() + ciBasis.getInternuclear_repulsion();
    //compare
    BOOST_CHECK(std::abs(en - (-8.0029560313)) < 1.0e-06); //energy from lih

}

/*
BOOST_AUTO_TEST_CASE ( doci_ref_cn_test_sparse ) {
    //cn_cation
    const std::string doci = "../tests/reference_data/doci_ref/cn_cation_ap1rog_631g.txt";
    //loading the basis
    doci::CI_basis ciBasis(doci);
    //calculating doci
    doci::SparseDOCI doci_test(ciBasis);
    //getting the electronic ground energy
    doci::State ground = doci_test.getGroundstates().at(0);
    //getting the total energy
    double en = ground.getEval() + ciBasis.getInternuclear_repulsion();
    //compare
    BOOST_CHECK(std::abs(en - (-91.7949190681)) < 1.0e-06); //energy from cn_cation

}

BOOST_AUTO_TEST_CASE ( doci_ref_co_test ) {
    //co_cation
    const std::string doci = "../tests/reference_data/doci_ref/co_cation_ap1rog_631g.txt";
    //loading the basis
    doci::CI_basis ciBasis(doci);
    //calculating doci
    doci::SparseDOCI doci_test(ciBasis);
    //getting the electronic ground energy
    doci::State ground = doci_test.getGroundstates().at(0);
    //getting the total energy
    double en = ground.getEval() + ciBasis.getInternuclear_repulsion();
    //compare
    BOOST_CHECK(std::abs(en - (-112.8392190666)) < 1.0e-06); //energy from co_cation

}

BOOST_AUTO_TEST_CASE ( doci_ref_h2o_test ) {
    //h2o
    const std::string doci = "../tests/reference_data/doci_ref/h2o_ap1rog_631g.txt";
    //loading the basis
    doci::CI_basis ciBasis(doci);
    //calculating doci
    doci::SparseDOCI doci_test(ciBasis);
    //getting the electronic ground energy
    doci::State ground = doci_test.getGroundstates().at(0);
    //getting the total energy
    double en = ground.getEval() + ciBasis.getInternuclear_repulsion();
    //compare
    BOOST_CHECK(std::abs(en - (-75.8179284026)) < 1.0e-06); //energy from h2o

}

BOOST_AUTO_TEST_CASE ( doci_ref_hf_test ) {
    //hf
    const std::string doci = "../tests/reference_data/doci_ref/hf_ap1rog_631g.txt";
    //loading the basis
    doci::CI_basis ciBasis(doci);
    //calculating doci
    doci::SparseDOCI doci_test(ciBasis);
    //getting the electronic ground energy
    doci::State ground = doci_test.getGroundstates().at(0);
    //getting the total energy
    double en = ground.getEval() + ciBasis.getInternuclear_repulsion();
    //compare
    BOOST_CHECK(std::abs(en - (-100.0704844678)) < 1.0e-06); //energy from hf

}


BOOST_AUTO_TEST_CASE ( doci_ref_lif_test ) {
    //h2o
    const std::string doci = "../tests/reference_data/doci_ref/lif_ap1rog_631g.txt";
    //loading the basis
    doci::CI_basis ciBasis(doci);
    //calculating doci
    doci::SparseDOCI doci_test(ciBasis);
    //getting the electronic ground energy
    doci::State ground = doci_test.getGroundstates().at(0);
    //getting the total energy
    double en = ground.getEval() + ciBasis.getInternuclear_repulsion();
    //compare
    BOOST_CHECK(std::abs(en - (-106.932928)) < 1.0e-06); //energy from lif

}

BOOST_AUTO_TEST_CASE ( doci_ref_nh3_test ) {
    //nh3
    const std::string doci = "../tests/reference_data/doci_ref/nh3_ap1rog_631g.txt";
    //loading the basis
    doci::CI_basis ciBasis(doci);
    //calculating doci
    doci::SparseDOCI doci_test(ciBasis);
    //getting the electronic ground energy
    doci::State ground = doci_test.getGroundstates().at(0);
    //getting the total energy
    double en = ground.getEval() + ciBasis.getInternuclear_repulsion();
    //compare
    BOOST_CHECK(std::abs(en - (-56.2703067338)) < 1.0e-06); //energy from nh3

}
*/