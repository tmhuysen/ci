# In this CMake file, we will find all required packages


# Find the boost package - needed for unittests
find_package(Boost REQUIRED)

# Find Armadillo for linear algebra operations
find_package(Armadillo REQUIRED)

# Find Eigen3
find_package(Eigen3 REQUIRED NO_MODULE)

# Find hf for HF calculations
find_package(hf 2.0.0 REQUIRED)

# Find libwrp for basis set and int calc
find_package(libwrp 2.1.1 REQUIRED)