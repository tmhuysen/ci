# In this CMake file, we will find all required packages


# Find the Boost package - needed for unittests
find_package(Boost REQUIRED)

# Find Eigen3
find_package(Eigen3 REQUIRED NO_MODULE)

# Find hf
find_package(hf 2.0.2 REQUIRED)

# Find the libint integral wrapper (also includes support for integral transformations)
find_package(libwint 2.2.2 REQUIRED)

# Find bmqc for bitset manipulations
find_package(bmqc 0.1.0 REQUIRED)

# We will use our custom FindSpectra.cmake-file, so that we can use find_package(Spectra)
find_package(Spectra REQUIRED)

# Find numopt
find_package(numopt 0.2.0 REQUIRED)
