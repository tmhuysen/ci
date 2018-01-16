# In this CMake file, we will include header files and link to libraries for a given test source

# for each test:
# ... add the boost headers ...
target_include_directories(${TEST_NAME} PUBLIC ${Boost_INCLUDE_DIRS})

# ... add this project's library ...
target_include_directories(${TEST_NAME} PRIVATE ${PROJECT_INCLUDE_FOLDER})
target_link_libraries(${TEST_NAME} PRIVATE ${LIBRARY_NAME})

# ... add Armadillo ...
target_include_directories(${TEST_NAME} PUBLIC ${ARMADILLO_INCUDE_DIRS})
target_link_libraries(${TEST_NAME} PUBLIC ${ARMADILLO_LIBRARIES})

# ... add Eigen3 ...
target_include_directories(${TEST_NAME} PUBLIC ${EIGEN3_INCLUDE_DIRS})
target_link_libraries(${TEST_NAME} PUBLIC ${EIGEN3_LIBRARIES})


# Include hf
target_include_directories(${TEST_NAME} PUBLIC ${hf_INCLUDE_DIRS})
target_link_libraries(${TEST_NAME} PUBLIC hf)

#include libwrp
target_include_directories(${TEST_NAME} PUBLIC ${libwrp_INCLUDE_DIRS})
target_link_libraries(${TEST_NAME} PUBLIC libwrp)
#set(ENV{MKLROOT} /opt/intel/compilers_and_libraries_2018.1.126/mac/mkl)
target_include_directories(${TEST_NAME} PRIVATE $ENV{MKLROOT}/include)
target_link_libraries(${TEST_NAME} PRIVATE mkl_intel_lp64 mkl_intel_thread mkl_core iomp5 pthread m dl)
#target_link_libraries(${TEST_NAME} PUBLIC ${MKLROOT}/lib/libmkl_intel_lp64.a ${MKLROOT}/lib/libmkl_sequential.a ${MKLROOT}/lib/libmkl_core.a -lpthread -lm -ldl)
#target_link_libraries(${TEST_NAME} PUBLIC ${MKLROOT}/lib/libmkl_intel_ilp64.a ${MKLROOT}/lib/libmkl_sequential.a ${MKLROOT}/lib/libmkl_core.a -lpthread -lm -ldl)

#target_link_libraries(${TEST_NAME} PUBLIC ${MKLROOT}/lib/libmkl_intel_ilp64.a ${MKLROOT}/lib/libmkl_sequential.a ${MKLROOT}/lib/libmkl_core.a -lpthread -lm -ldl)