# doci
[![Build Status](https://travis-ci.org/GQCG/doci.svg?branch=master)](https://travis-ci.org/GQCG/doci)

A C++ library for performing DOCI calculations.

## Dependencies
[![Eigen3 Dependency](https://img.shields.io/badge/Eigen-3+-blue.svg)](http://eigen.tuxfamily.org/index.php?title=Main_Page)
[![libwint Dependency](https://img.shields.io/badge/libwrp-2.2.2+-blue.svg)](https://github.com/GQCG/libwrp)
[![hf Dependency](https://img.shields.io/badge/hf-2.0.2+-blue.svg)](https://github.com/GQCG/hf)

## Installation
To install this library:
1. clone the master branch

        https://github.com/GQCG/doci.git
        cd doci

2. perform an out-of-source cmake build:

        mkdir build && cd build
        cmake -DINSTALLATION_PREFIX=prefix ..
        make && make test && sudo make install

    where
    * `prefix` is the installation prefix (defaulted to `/usr/local`) you want the library to be installed at:
        * the library `libdoci.a` will be installed in `prefix/doci/lib`
        * the header files (and cmake files, see Usage) will be installed in `prefix/doci/include`


## Usage
Basic usage of this library can be found in the `tests` directory. If you use CMake in other projects, you can add the following CMake command to the CMakeLists.txt-file:

    find_package(doci x.y.z)

where `x.y.z` is the version number. CMake then provides the commands `doci_INCLUDE_DIRS` to be used in your `target_include_directories` and the library `doci` to be used in your `target_link_libraries`.
