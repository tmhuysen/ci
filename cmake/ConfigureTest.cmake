# In this CMake file, we will include header files and link to libraries for a given test source

# for each test:
# ... add the boost headers ...
target_include_directories(${TEST_NAME} PUBLIC ${Boost_INCLUDE_DIRS})

# ... add this project's library ...
target_include_directories(${TEST_NAME} PRIVATE ${PROJECT_INCLUDE_FOLDER})
target_link_libraries(${TEST_NAME} PRIVATE ${LIBRARY_NAME})

# ... add numopt ...
target_include_directories(${TEST_NAME} PUBLIC ${numopt_INCLUDE_DIRS})
target_link_libraries(${TEST_NAME} PUBLIC numopt)

# ... add Eigen3 ...
target_link_libraries(${TEST_NAME} PUBLIC Eigen3::Eigen)

# ... add Spectra ...
target_include_directories(${TEST_NAME} PUBLIC /opt/local/spectra/include)

# ... include hf ...
target_include_directories(${TEST_NAME} PUBLIC ${hf_INCLUDE_DIRS})
target_link_libraries(${TEST_NAME} PUBLIC hf)

# ... add bmqc ...
target_include_directories(${TEST_NAME} PUBLIC ${bmqc_INCLUDE_DIRS})
target_link_libraries(${TEST_NAME} PUBLIC bmqc)


