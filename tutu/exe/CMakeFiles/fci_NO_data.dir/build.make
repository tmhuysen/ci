# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.11

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
CMAKE_SOURCE_DIR = /Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/Forks/ci

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/Forks/ci/tutu

# Include any dependencies generated for this target.
include exe/CMakeFiles/fci_NO_data.dir/depend.make

# Include the progress variables for this target.
include exe/CMakeFiles/fci_NO_data.dir/progress.make

# Include the compile flags for this target's objects.
include exe/CMakeFiles/fci_NO_data.dir/flags.make

exe/CMakeFiles/fci_NO_data.dir/fci_NO_data.cpp.o: exe/CMakeFiles/fci_NO_data.dir/flags.make
exe/CMakeFiles/fci_NO_data.dir/fci_NO_data.cpp.o: ../exe/fci_NO_data.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/Forks/ci/tutu/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object exe/CMakeFiles/fci_NO_data.dir/fci_NO_data.cpp.o"
	cd /Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/Forks/ci/tutu/exe && /opt/local/bin/clang++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/fci_NO_data.dir/fci_NO_data.cpp.o -c /Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/Forks/ci/exe/fci_NO_data.cpp

exe/CMakeFiles/fci_NO_data.dir/fci_NO_data.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/fci_NO_data.dir/fci_NO_data.cpp.i"
	cd /Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/Forks/ci/tutu/exe && /opt/local/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/Forks/ci/exe/fci_NO_data.cpp > CMakeFiles/fci_NO_data.dir/fci_NO_data.cpp.i

exe/CMakeFiles/fci_NO_data.dir/fci_NO_data.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/fci_NO_data.dir/fci_NO_data.cpp.s"
	cd /Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/Forks/ci/tutu/exe && /opt/local/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/Forks/ci/exe/fci_NO_data.cpp -o CMakeFiles/fci_NO_data.dir/fci_NO_data.cpp.s

# Object files for target fci_NO_data
fci_NO_data_OBJECTS = \
"CMakeFiles/fci_NO_data.dir/fci_NO_data.cpp.o"

# External object files for target fci_NO_data
fci_NO_data_EXTERNAL_OBJECTS =

exe/fci_NO_data: exe/CMakeFiles/fci_NO_data.dir/fci_NO_data.cpp.o
exe/fci_NO_data: exe/CMakeFiles/fci_NO_data.dir/build.make
exe/fci_NO_data: src/libci.a
exe/fci_NO_data: /usr/local/hf/lib/libhf.a
exe/fci_NO_data: /usr/local/bmqc/lib/libbmqc.a
exe/fci_NO_data: /usr/local/cpputil/lib/libcpputil.a
exe/fci_NO_data: /usr/local/libwint/lib/libwint.a
exe/fci_NO_data: /usr/local/libint/2.3.1/lib/libint2.a
exe/fci_NO_data: /usr/local/numopt/lib/libnumopt.a
exe/fci_NO_data: /opt/local/lib/libarmadillo.dylib
exe/fci_NO_data: exe/CMakeFiles/fci_NO_data.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/Forks/ci/tutu/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable fci_NO_data"
	cd /Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/Forks/ci/tutu/exe && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/fci_NO_data.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
exe/CMakeFiles/fci_NO_data.dir/build: exe/fci_NO_data

.PHONY : exe/CMakeFiles/fci_NO_data.dir/build

exe/CMakeFiles/fci_NO_data.dir/clean:
	cd /Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/Forks/ci/tutu/exe && $(CMAKE_COMMAND) -P CMakeFiles/fci_NO_data.dir/cmake_clean.cmake
.PHONY : exe/CMakeFiles/fci_NO_data.dir/clean

exe/CMakeFiles/fci_NO_data.dir/depend:
	cd /Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/Forks/ci/tutu && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/Forks/ci /Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/Forks/ci/exe /Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/Forks/ci/tutu /Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/Forks/ci/tutu/exe /Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/Forks/ci/tutu/exe/CMakeFiles/fci_NO_data.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : exe/CMakeFiles/fci_NO_data.dir/depend

