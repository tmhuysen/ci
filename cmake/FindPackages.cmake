# In this CMake file, we will find all required packages


# Find the Boost package - needed for unittests
find_package(Boost REQUIRED)

# Find Eigen3
find_package(Eigen3 REQUIRED NO_MODULE)

# Find hf
find_package(hf 2.0.2 REQUIRED)

# Find the libint integral wrapper + support for integral transformations
find_package(libwint 2.2.2 REQUIRED)
