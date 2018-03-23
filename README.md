# ci v1.0.0
[![Build Status](https://travis-ci.org/GQCG/ci.svg?branch=master)](https://travis-ci.org/GQCG/ci)

A C++ library for performing configuration interaction (CI) calculations.

## Dependencies
[![Eigen3 Dependency](https://img.shields.io/badge/Eigen-3+-blue.svg)](http://eigen.tuxfamily.org/index.php?title=Main_Page)
[![bmqc Dependency](https://img.shields.io/badge/bmqc-1.0.1+-blue.svg)](https://github.com/GQCG/bmqc)
[![libwint Dependency](https://img.shields.io/badge/libwrp-3.0.0+-blue.svg)](https://github.com/GQCG/libwrp)
[![hf Dependency](https://img.shields.io/badge/hf-3.0.0+-blue.svg)](https://github.com/GQCG/hf)
[![numopt Dependency](https://img.shields.io/badge/bmqc-1.0.0+-blue.svg)](https://github.com/GQCG/numopt)


## Installation
To install this library:
1. clone the master branch, which contains the latest release

        https://github.com/GQCG/ci.git --branch master --single-branch
        cd ci

2. perform an out-of-source cmake build:

        mkdir build && cd build
        cmake -DINSTALLATION_PREFIX=prefix ..
        make && make test && sudo make install

    where
    * `prefix` is the installation prefix (defaulted to `/usr/local`) you want the library to be installed at:
        * the library `libci.a` will be installed in `prefix/ci/lib`
        * the header files (and cmake files, see Usage) will be installed in `prefix/ci/include`


## Usage
Basic usage of this library can be found in the `tests` directory. If you use CMake in other projects, you can add the following CMake command to the CMakeLists.txt-file:

    find_package(ci x.y.z)

where `x.y.z` is the version number. CMake then provides the commands `ci_INCLUDE_DIRS` to be used in your `target_include_directories` and the library `ci` to be used in your `target_link_libraries`.
