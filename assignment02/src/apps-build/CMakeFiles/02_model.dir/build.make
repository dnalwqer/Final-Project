# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.3

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
CMAKE_COMMAND = /Applications/CMake.app/Contents/bin/cmake

# The command to remove a file.
RM = /Applications/CMake.app/Contents/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/dnalwqer/Desktop/CG/assignment02/src/apps

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/dnalwqer/Desktop/CG/assignment02/src/apps-build

# Include any dependencies generated for this target.
include CMakeFiles/02_model.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/02_model.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/02_model.dir/flags.make

CMakeFiles/02_model.dir/02_model.o: CMakeFiles/02_model.dir/flags.make
CMakeFiles/02_model.dir/02_model.o: /Users/dnalwqer/Desktop/CG/assignment02/src/apps/02_model.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/dnalwqer/Desktop/CG/assignment02/src/apps-build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/02_model.dir/02_model.o"
	/usr/bin/clang++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/02_model.dir/02_model.o -c /Users/dnalwqer/Desktop/CG/assignment02/src/apps/02_model.cpp

CMakeFiles/02_model.dir/02_model.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/02_model.dir/02_model.i"
	/usr/bin/clang++  $(CXX_DEFINES) $(CXX_FLAGS) -E /Users/dnalwqer/Desktop/CG/assignment02/src/apps/02_model.cpp > CMakeFiles/02_model.dir/02_model.i

CMakeFiles/02_model.dir/02_model.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/02_model.dir/02_model.s"
	/usr/bin/clang++  $(CXX_DEFINES) $(CXX_FLAGS) -S /Users/dnalwqer/Desktop/CG/assignment02/src/apps/02_model.cpp -o CMakeFiles/02_model.dir/02_model.s

CMakeFiles/02_model.dir/02_model.o.requires:

.PHONY : CMakeFiles/02_model.dir/02_model.o.requires

CMakeFiles/02_model.dir/02_model.o.provides: CMakeFiles/02_model.dir/02_model.o.requires
	$(MAKE) -f CMakeFiles/02_model.dir/build.make CMakeFiles/02_model.dir/02_model.o.provides.build
.PHONY : CMakeFiles/02_model.dir/02_model.o.provides

CMakeFiles/02_model.dir/02_model.o.provides.build: CMakeFiles/02_model.dir/02_model.o


# Object files for target 02_model
02_model_OBJECTS = \
"CMakeFiles/02_model.dir/02_model.o"

# External object files for target 02_model
02_model_EXTERNAL_OBJECTS =

02_model: CMakeFiles/02_model.dir/02_model.o
02_model: CMakeFiles/02_model.dir/build.make
02_model: CMakeFiles/02_model.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/dnalwqer/Desktop/CG/assignment02/src/apps-build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable 02_model"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/02_model.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/02_model.dir/build: 02_model

.PHONY : CMakeFiles/02_model.dir/build

CMakeFiles/02_model.dir/requires: CMakeFiles/02_model.dir/02_model.o.requires

.PHONY : CMakeFiles/02_model.dir/requires

CMakeFiles/02_model.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/02_model.dir/cmake_clean.cmake
.PHONY : CMakeFiles/02_model.dir/clean

CMakeFiles/02_model.dir/depend:
	cd /Users/dnalwqer/Desktop/CG/assignment02/src/apps-build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/dnalwqer/Desktop/CG/assignment02/src/apps /Users/dnalwqer/Desktop/CG/assignment02/src/apps /Users/dnalwqer/Desktop/CG/assignment02/src/apps-build /Users/dnalwqer/Desktop/CG/assignment02/src/apps-build /Users/dnalwqer/Desktop/CG/assignment02/src/apps-build/CMakeFiles/02_model.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/02_model.dir/depend
