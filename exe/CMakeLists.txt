# exe-level CMakeLists

# Add an executable for davidson timings
add_executable(fci_NO_data ${PROJECT_EXECUTABLES_FOLDER}/fci_NO_data.cpp)

# Add this project's headers
configure_executable(fci_NO_data)
