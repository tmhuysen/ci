# In this CMake file, we will find all required packages


# Find the Boost package - needed for unittests
find_package(Boost REQUIRED)

# Find numopt for DavidsonSolver
find_package(numopt 0.2.0 REQUIRED)

# Find Eigen3
find_package(Eigen3 REQUIRED NO_MODULE)

# Find hf
find_package(hf 2.0.2 REQUIRED)

# Find the libint integral wrapper (also includes support for integral transformations)
find_package(libwint 2.2.2 REQUIRED)

# Find bmqc for bitset manipulations
find_package(bmqc 0.1.0 REQUIRED)


# Find numopt for DavidsonSolver
find_package(numopt 0.2.0 REQUIRED)