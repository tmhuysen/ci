# src-level CMakeLists.txt
cmake_minimum_required(VERSION 3.8)


project(ci
        VERSION 0.0.1
        LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)



# Install FindSpectra.cmake. This is done before FindPackages.cmake, since we are using find_package(Spectra)
# CMake looks for Find<package>.cmake files in its ${CMAKE_MODULE_PATH}. By default, this is empty. Since we don't want to bother with installing our custom-made FindSpectra.cmake, we will set the ${CMAKE_MODULE_PATH} accordingly
set(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH};${CMAKE_SOURCE_DIR}/cmake")
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

# Set CMake variables
include(${CMAKE_SOURCE_DIR}/cmake/SetCMakeEnvironment.cmake)

# Find all required packages
include(${CMAKE_SOURCE_DIR}/cmake/FindPackages.cmake)

# Add the sources for the library
add_subdirectory(${PROJECT_SOURCE_FOLDER})

# Set up testing
enable_testing()
add_subdirectory(${PROJECT_TESTS_FOLDER})
