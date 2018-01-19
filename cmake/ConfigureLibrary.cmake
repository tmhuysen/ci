# In this CMake file, we will include the headers and link to the necessary libraries


# Include this project's headers
target_include_directories(${LIBRARY_NAME} PRIVATE ${PROJECT_INCLUDE_FOLDER})

# Include the boost headers
target_include_directories(${LIBRARY_NAME} PUBLIC ${Boost_INCLUDE_DIRS})

# Include Armadillo
target_include_directories(${LIBRARY_NAME} PUBLIC ${ARMADILLO_INCUDE_DIRS})
target_link_libraries(${LIBRARY_NAME} PUBLIC ${ARMADILLO_LIBRARIES})

# Include Eigen
target_include_directories(${LIBRARY_NAME} PUBLIC ${Eigen3_INCUDE_DIRS})
target_link_libraries(${LIBRARY_NAME} PUBLIC Eigen3::Eigen)

# Include hf
target_include_directories(${LIBRARY_NAME} PUBLIC ${hf_INCLUDE_DIRS})
target_link_libraries(${LIBRARY_NAME} PUBLIC hf)

#include libwrp
target_include_directories(${LIBRARY_NAME} PUBLIC ${libwrp_INCLUDE_DIRS})
target_link_libraries(${LIBRARY_NAME} PUBLIC libwrp)
#set(ENV{MKLROOT} /opt/intel/compilers_and_libraries_2018.1.126/mac/mkl)
target_include_directories(${LIBRARY_NAME} PRIVATE $ENV{MKLROOT}/include)
#target_link_libraries(${LIBRARY_NAME} PRIVATE mkl_intel_lp64 mkl_intel_thread mkl_core iomp5 pthread m dl)
target_link_libraries(${LIBRARY_NAME} PRIVATE $ENV{MKLROOT}/lib/libmkl_intel_lp64.a $ENV{MKLROOT}/lib/libmkl_sequential.a $ENV{MKLROOT}/lib/libmkl_core.a -lpthread -lm -ldl)
#target_link_libraries(${LIBRARY_NAME} PUBLIC ${MKLROOT}/lib/libmkl_intel_ilp64.a ${MKLROOT}/lib/libmkl_sequential.a ${MKLROOT}/lib/libmkl_core.a -lpthread -lm -ldl)
# -m64 -I${MKLROOT}/include
#  -DMKL_ILP64 -m64 -I${MKLROOT}/include