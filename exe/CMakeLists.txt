# exe-level CMakeLists

# Add an executable for davidson timings
add_executable(fci_NO_data ${PROJECT_EXECUTABLES_FOLDER}/fci_NO_data.cpp)
add_executable(fci_NO_data3 ${PROJECT_EXECUTABLES_FOLDER}/fci_NO_data3.cpp)
add_executable(10Bohr ${PROJECT_EXECUTABLES_FOLDER}/10Bohr.cpp)

# Add this project's headers
configure_executable(fci_NO_data)
configure_executable(fci_NO_data3)
configure_executable(10Bohr)
