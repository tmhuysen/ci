# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.10

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /opt/local/bin/cmake

# The command to remove a file.
RM = /opt/local/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/GQCG/ci

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/GQCG/ci/legit_test

# Include any dependencies generated for this target.
include tests/CMakeFiles/Hamiltonian_test.dir/depend.make

# Include the progress variables for this target.
include tests/CMakeFiles/Hamiltonian_test.dir/progress.make

# Include the compile flags for this target's objects.
include tests/CMakeFiles/Hamiltonian_test.dir/flags.make

tests/CMakeFiles/Hamiltonian_test.dir/Hamiltonian_test.cpp.o: tests/CMakeFiles/Hamiltonian_test.dir/flags.make
tests/CMakeFiles/Hamiltonian_test.dir/Hamiltonian_test.cpp.o: ../tests/Hamiltonian_test.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/GQCG/ci/legit_test/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object tests/CMakeFiles/Hamiltonian_test.dir/Hamiltonian_test.cpp.o"
	cd /Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/GQCG/ci/legit_test/tests && /opt/local/bin/clang++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/Hamiltonian_test.dir/Hamiltonian_test.cpp.o -c /Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/GQCG/ci/tests/Hamiltonian_test.cpp

tests/CMakeFiles/Hamiltonian_test.dir/Hamiltonian_test.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/Hamiltonian_test.dir/Hamiltonian_test.cpp.i"
	cd /Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/GQCG/ci/legit_test/tests && /opt/local/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/GQCG/ci/tests/Hamiltonian_test.cpp > CMakeFiles/Hamiltonian_test.dir/Hamiltonian_test.cpp.i

tests/CMakeFiles/Hamiltonian_test.dir/Hamiltonian_test.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/Hamiltonian_test.dir/Hamiltonian_test.cpp.s"
	cd /Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/GQCG/ci/legit_test/tests && /opt/local/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/GQCG/ci/tests/Hamiltonian_test.cpp -o CMakeFiles/Hamiltonian_test.dir/Hamiltonian_test.cpp.s

tests/CMakeFiles/Hamiltonian_test.dir/Hamiltonian_test.cpp.o.requires:

.PHONY : tests/CMakeFiles/Hamiltonian_test.dir/Hamiltonian_test.cpp.o.requires

tests/CMakeFiles/Hamiltonian_test.dir/Hamiltonian_test.cpp.o.provides: tests/CMakeFiles/Hamiltonian_test.dir/Hamiltonian_test.cpp.o.requires
	$(MAKE) -f tests/CMakeFiles/Hamiltonian_test.dir/build.make tests/CMakeFiles/Hamiltonian_test.dir/Hamiltonian_test.cpp.o.provides.build
.PHONY : tests/CMakeFiles/Hamiltonian_test.dir/Hamiltonian_test.cpp.o.provides

tests/CMakeFiles/Hamiltonian_test.dir/Hamiltonian_test.cpp.o.provides.build: tests/CMakeFiles/Hamiltonian_test.dir/Hamiltonian_test.cpp.o


# Object files for target Hamiltonian_test
Hamiltonian_test_OBJECTS = \
"CMakeFiles/Hamiltonian_test.dir/Hamiltonian_test.cpp.o"

# External object files for target Hamiltonian_test
Hamiltonian_test_EXTERNAL_OBJECTS =

tests/Hamiltonian_test: tests/CMakeFiles/Hamiltonian_test.dir/Hamiltonian_test.cpp.o
tests/Hamiltonian_test: tests/CMakeFiles/Hamiltonian_test.dir/build.make
tests/Hamiltonian_test: src/libci.a
tests/Hamiltonian_test: /usr/local/numopt/lib/libnumopt.a
tests/Hamiltonian_test: /usr/local/hf/lib/libhf.a
tests/Hamiltonian_test: /usr/local/bmqc/lib/libbmqc.a
tests/Hamiltonian_test: /opt/intel/compilers_and_libraries_2018.1.126/mac/mkl/lib/libmkl_intel_lp64.a
tests/Hamiltonian_test: /opt/intel/compilers_and_libraries_2018.1.126/mac/mkl/lib/libmkl_sequential.a
tests/Hamiltonian_test: /opt/intel/compilers_and_libraries_2018.1.126/mac/mkl/lib/libmkl_core.a
tests/Hamiltonian_test: /usr/local/libwint/lib/libwint.a
tests/Hamiltonian_test: /usr/local/libint/2.3.1/lib/libint2.a
tests/Hamiltonian_test: tests/CMakeFiles/Hamiltonian_test.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/GQCG/ci/legit_test/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable Hamiltonian_test"
	cd /Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/GQCG/ci/legit_test/tests && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/Hamiltonian_test.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
tests/CMakeFiles/Hamiltonian_test.dir/build: tests/Hamiltonian_test

.PHONY : tests/CMakeFiles/Hamiltonian_test.dir/build

tests/CMakeFiles/Hamiltonian_test.dir/requires: tests/CMakeFiles/Hamiltonian_test.dir/Hamiltonian_test.cpp.o.requires

.PHONY : tests/CMakeFiles/Hamiltonian_test.dir/requires

tests/CMakeFiles/Hamiltonian_test.dir/clean:
	cd /Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/GQCG/ci/legit_test/tests && $(CMAKE_COMMAND) -P CMakeFiles/Hamiltonian_test.dir/cmake_clean.cmake
.PHONY : tests/CMakeFiles/Hamiltonian_test.dir/clean

tests/CMakeFiles/Hamiltonian_test.dir/depend:
	cd /Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/GQCG/ci/legit_test && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/GQCG/ci /Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/GQCG/ci/tests /Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/GQCG/ci/legit_test /Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/GQCG/ci/legit_test/tests /Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/GQCG/ci/legit_test/tests/CMakeFiles/Hamiltonian_test.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : tests/CMakeFiles/Hamiltonian_test.dir/depend

