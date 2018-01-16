# Configuration file for the "DOCI_Head" package
# It defines the following variables:
#   doci_head_INCLUDE_DIRS   - header include directories
#   doci_head_LIBRARIES      - library to link against


# Specify the include directory
set(doci_head_INCLUDE_DIRS /usr/local/doci_head/include)

# Import the exported targets
include(/usr/local/doci_head/cmake/DOCI_HeadTargets.cmake)

# Specify the library value
set(doci_head_LIBRARIES DOCI_Head)

# Output the library version
message(STATUS "doci_head version: 0.0.1")
