# In this CMake file, we will include header files and link to libraries for a given test source

# for each test:
# ... add the boost headers ...
target_include_directories(${TEST_NAME} PUBLIC ${Boost_INCLUDE_DIRS})

# ... add this project's library ...
target_include_directories(${TEST_NAME} PUBLIC ${PROJECT_INCLUDE_FOLDER})
target_link_libraries(${TEST_NAME} PUBLIC ${LIBRARY_NAME})

# ... include hf ...
target_include_directories(${TEST_NAME} PUBLIC ${hf_INCLUDE_DIRS})
target_link_libraries(${TEST_NAME} PUBLIC hf)
