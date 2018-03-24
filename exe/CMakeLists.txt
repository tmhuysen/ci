# exe-level CMakeLists

# Add an executable for davidson timings
add_executable(dense_sparse_DOCI_timings ${PROJECT_EXECUTABLES_FOLDER}/dense_sparse_DOCI_timings.cpp)

# Add this project's headers
configure_executable(dense_sparse_DOCI_timings)
