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
CMAKE_SOURCE_DIR = /Users/dnalwqer/Desktop/CG/Dart-Graphics-Assign01-master

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/dnalwqer/Desktop/CG/Dart-Graphics-Assign01-master-build

# Include any dependencies generated for this target.
include src/common/CMakeFiles/common.dir/depend.make

# Include the progress variables for this target.
include src/common/CMakeFiles/common.dir/progress.make

# Include the compile flags for this target's objects.
include src/common/CMakeFiles/common.dir/flags.make

src/common/CMakeFiles/common.dir/image.cpp.o: src/common/CMakeFiles/common.dir/flags.make
src/common/CMakeFiles/common.dir/image.cpp.o: /Users/dnalwqer/Desktop/CG/Dart-Graphics-Assign01-master/src/common/image.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/dnalwqer/Desktop/CG/Dart-Graphics-Assign01-master-build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object src/common/CMakeFiles/common.dir/image.cpp.o"
	cd /Users/dnalwqer/Desktop/CG/Dart-Graphics-Assign01-master-build/src/common && /usr/bin/clang++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/common.dir/image.cpp.o -c /Users/dnalwqer/Desktop/CG/Dart-Graphics-Assign01-master/src/common/image.cpp

src/common/CMakeFiles/common.dir/image.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/common.dir/image.cpp.i"
	cd /Users/dnalwqer/Desktop/CG/Dart-Graphics-Assign01-master-build/src/common && /usr/bin/clang++  $(CXX_DEFINES) $(CXX_FLAGS) -E /Users/dnalwqer/Desktop/CG/Dart-Graphics-Assign01-master/src/common/image.cpp > CMakeFiles/common.dir/image.cpp.i

src/common/CMakeFiles/common.dir/image.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/common.dir/image.cpp.s"
	cd /Users/dnalwqer/Desktop/CG/Dart-Graphics-Assign01-master-build/src/common && /usr/bin/clang++  $(CXX_DEFINES) $(CXX_FLAGS) -S /Users/dnalwqer/Desktop/CG/Dart-Graphics-Assign01-master/src/common/image.cpp -o CMakeFiles/common.dir/image.cpp.s

src/common/CMakeFiles/common.dir/image.cpp.o.requires:

.PHONY : src/common/CMakeFiles/common.dir/image.cpp.o.requires

src/common/CMakeFiles/common.dir/image.cpp.o.provides: src/common/CMakeFiles/common.dir/image.cpp.o.requires
	$(MAKE) -f src/common/CMakeFiles/common.dir/build.make src/common/CMakeFiles/common.dir/image.cpp.o.provides.build
.PHONY : src/common/CMakeFiles/common.dir/image.cpp.o.provides

src/common/CMakeFiles/common.dir/image.cpp.o.provides.build: src/common/CMakeFiles/common.dir/image.cpp.o


src/common/CMakeFiles/common.dir/json.cpp.o: src/common/CMakeFiles/common.dir/flags.make
src/common/CMakeFiles/common.dir/json.cpp.o: /Users/dnalwqer/Desktop/CG/Dart-Graphics-Assign01-master/src/common/json.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/dnalwqer/Desktop/CG/Dart-Graphics-Assign01-master-build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object src/common/CMakeFiles/common.dir/json.cpp.o"
	cd /Users/dnalwqer/Desktop/CG/Dart-Graphics-Assign01-master-build/src/common && /usr/bin/clang++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/common.dir/json.cpp.o -c /Users/dnalwqer/Desktop/CG/Dart-Graphics-Assign01-master/src/common/json.cpp

src/common/CMakeFiles/common.dir/json.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/common.dir/json.cpp.i"
	cd /Users/dnalwqer/Desktop/CG/Dart-Graphics-Assign01-master-build/src/common && /usr/bin/clang++  $(CXX_DEFINES) $(CXX_FLAGS) -E /Users/dnalwqer/Desktop/CG/Dart-Graphics-Assign01-master/src/common/json.cpp > CMakeFiles/common.dir/json.cpp.i

src/common/CMakeFiles/common.dir/json.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/common.dir/json.cpp.s"
	cd /Users/dnalwqer/Desktop/CG/Dart-Graphics-Assign01-master-build/src/common && /usr/bin/clang++  $(CXX_DEFINES) $(CXX_FLAGS) -S /Users/dnalwqer/Desktop/CG/Dart-Graphics-Assign01-master/src/common/json.cpp -o CMakeFiles/common.dir/json.cpp.s

src/common/CMakeFiles/common.dir/json.cpp.o.requires:

.PHONY : src/common/CMakeFiles/common.dir/json.cpp.o.requires

src/common/CMakeFiles/common.dir/json.cpp.o.provides: src/common/CMakeFiles/common.dir/json.cpp.o.requires
	$(MAKE) -f src/common/CMakeFiles/common.dir/build.make src/common/CMakeFiles/common.dir/json.cpp.o.provides.build
.PHONY : src/common/CMakeFiles/common.dir/json.cpp.o.provides

src/common/CMakeFiles/common.dir/json.cpp.o.provides.build: src/common/CMakeFiles/common.dir/json.cpp.o


src/common/CMakeFiles/common.dir/scene.cpp.o: src/common/CMakeFiles/common.dir/flags.make
src/common/CMakeFiles/common.dir/scene.cpp.o: /Users/dnalwqer/Desktop/CG/Dart-Graphics-Assign01-master/src/common/scene.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/dnalwqer/Desktop/CG/Dart-Graphics-Assign01-master-build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object src/common/CMakeFiles/common.dir/scene.cpp.o"
	cd /Users/dnalwqer/Desktop/CG/Dart-Graphics-Assign01-master-build/src/common && /usr/bin/clang++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/common.dir/scene.cpp.o -c /Users/dnalwqer/Desktop/CG/Dart-Graphics-Assign01-master/src/common/scene.cpp

src/common/CMakeFiles/common.dir/scene.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/common.dir/scene.cpp.i"
	cd /Users/dnalwqer/Desktop/CG/Dart-Graphics-Assign01-master-build/src/common && /usr/bin/clang++  $(CXX_DEFINES) $(CXX_FLAGS) -E /Users/dnalwqer/Desktop/CG/Dart-Graphics-Assign01-master/src/common/scene.cpp > CMakeFiles/common.dir/scene.cpp.i

src/common/CMakeFiles/common.dir/scene.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/common.dir/scene.cpp.s"
	cd /Users/dnalwqer/Desktop/CG/Dart-Graphics-Assign01-master-build/src/common && /usr/bin/clang++  $(CXX_DEFINES) $(CXX_FLAGS) -S /Users/dnalwqer/Desktop/CG/Dart-Graphics-Assign01-master/src/common/scene.cpp -o CMakeFiles/common.dir/scene.cpp.s

src/common/CMakeFiles/common.dir/scene.cpp.o.requires:

.PHONY : src/common/CMakeFiles/common.dir/scene.cpp.o.requires

src/common/CMakeFiles/common.dir/scene.cpp.o.provides: src/common/CMakeFiles/common.dir/scene.cpp.o.requires
	$(MAKE) -f src/common/CMakeFiles/common.dir/build.make src/common/CMakeFiles/common.dir/scene.cpp.o.provides.build
.PHONY : src/common/CMakeFiles/common.dir/scene.cpp.o.provides

src/common/CMakeFiles/common.dir/scene.cpp.o.provides.build: src/common/CMakeFiles/common.dir/scene.cpp.o


src/common/CMakeFiles/common.dir/ext/lodepng/lodepng.cpp.o: src/common/CMakeFiles/common.dir/flags.make
src/common/CMakeFiles/common.dir/ext/lodepng/lodepng.cpp.o: /Users/dnalwqer/Desktop/CG/Dart-Graphics-Assign01-master/src/common/ext/lodepng/lodepng.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/dnalwqer/Desktop/CG/Dart-Graphics-Assign01-master-build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object src/common/CMakeFiles/common.dir/ext/lodepng/lodepng.cpp.o"
	cd /Users/dnalwqer/Desktop/CG/Dart-Graphics-Assign01-master-build/src/common && /usr/bin/clang++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/common.dir/ext/lodepng/lodepng.cpp.o -c /Users/dnalwqer/Desktop/CG/Dart-Graphics-Assign01-master/src/common/ext/lodepng/lodepng.cpp

src/common/CMakeFiles/common.dir/ext/lodepng/lodepng.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/common.dir/ext/lodepng/lodepng.cpp.i"
	cd /Users/dnalwqer/Desktop/CG/Dart-Graphics-Assign01-master-build/src/common && /usr/bin/clang++  $(CXX_DEFINES) $(CXX_FLAGS) -E /Users/dnalwqer/Desktop/CG/Dart-Graphics-Assign01-master/src/common/ext/lodepng/lodepng.cpp > CMakeFiles/common.dir/ext/lodepng/lodepng.cpp.i

src/common/CMakeFiles/common.dir/ext/lodepng/lodepng.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/common.dir/ext/lodepng/lodepng.cpp.s"
	cd /Users/dnalwqer/Desktop/CG/Dart-Graphics-Assign01-master-build/src/common && /usr/bin/clang++  $(CXX_DEFINES) $(CXX_FLAGS) -S /Users/dnalwqer/Desktop/CG/Dart-Graphics-Assign01-master/src/common/ext/lodepng/lodepng.cpp -o CMakeFiles/common.dir/ext/lodepng/lodepng.cpp.s

src/common/CMakeFiles/common.dir/ext/lodepng/lodepng.cpp.o.requires:

.PHONY : src/common/CMakeFiles/common.dir/ext/lodepng/lodepng.cpp.o.requires

src/common/CMakeFiles/common.dir/ext/lodepng/lodepng.cpp.o.provides: src/common/CMakeFiles/common.dir/ext/lodepng/lodepng.cpp.o.requires
	$(MAKE) -f src/common/CMakeFiles/common.dir/build.make src/common/CMakeFiles/common.dir/ext/lodepng/lodepng.cpp.o.provides.build
.PHONY : src/common/CMakeFiles/common.dir/ext/lodepng/lodepng.cpp.o.provides

src/common/CMakeFiles/common.dir/ext/lodepng/lodepng.cpp.o.provides.build: src/common/CMakeFiles/common.dir/ext/lodepng/lodepng.cpp.o


src/common/CMakeFiles/common.dir/ext/glew/glew.c.o: src/common/CMakeFiles/common.dir/flags.make
src/common/CMakeFiles/common.dir/ext/glew/glew.c.o: /Users/dnalwqer/Desktop/CG/Dart-Graphics-Assign01-master/src/common/ext/glew/glew.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/dnalwqer/Desktop/CG/Dart-Graphics-Assign01-master-build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building C object src/common/CMakeFiles/common.dir/ext/glew/glew.c.o"
	cd /Users/dnalwqer/Desktop/CG/Dart-Graphics-Assign01-master-build/src/common && /usr/bin/clang  $(C_DEFINES) $(C_FLAGS) -o CMakeFiles/common.dir/ext/glew/glew.c.o   -c /Users/dnalwqer/Desktop/CG/Dart-Graphics-Assign01-master/src/common/ext/glew/glew.c

src/common/CMakeFiles/common.dir/ext/glew/glew.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/common.dir/ext/glew/glew.c.i"
	cd /Users/dnalwqer/Desktop/CG/Dart-Graphics-Assign01-master-build/src/common && /usr/bin/clang  $(C_DEFINES) $(C_FLAGS) -E /Users/dnalwqer/Desktop/CG/Dart-Graphics-Assign01-master/src/common/ext/glew/glew.c > CMakeFiles/common.dir/ext/glew/glew.c.i

src/common/CMakeFiles/common.dir/ext/glew/glew.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/common.dir/ext/glew/glew.c.s"
	cd /Users/dnalwqer/Desktop/CG/Dart-Graphics-Assign01-master-build/src/common && /usr/bin/clang  $(C_DEFINES) $(C_FLAGS) -S /Users/dnalwqer/Desktop/CG/Dart-Graphics-Assign01-master/src/common/ext/glew/glew.c -o CMakeFiles/common.dir/ext/glew/glew.c.s

src/common/CMakeFiles/common.dir/ext/glew/glew.c.o.requires:

.PHONY : src/common/CMakeFiles/common.dir/ext/glew/glew.c.o.requires

src/common/CMakeFiles/common.dir/ext/glew/glew.c.o.provides: src/common/CMakeFiles/common.dir/ext/glew/glew.c.o.requires
	$(MAKE) -f src/common/CMakeFiles/common.dir/build.make src/common/CMakeFiles/common.dir/ext/glew/glew.c.o.provides.build
.PHONY : src/common/CMakeFiles/common.dir/ext/glew/glew.c.o.provides

src/common/CMakeFiles/common.dir/ext/glew/glew.c.o.provides.build: src/common/CMakeFiles/common.dir/ext/glew/glew.c.o


# Object files for target common
common_OBJECTS = \
"CMakeFiles/common.dir/image.cpp.o" \
"CMakeFiles/common.dir/json.cpp.o" \
"CMakeFiles/common.dir/scene.cpp.o" \
"CMakeFiles/common.dir/ext/lodepng/lodepng.cpp.o" \
"CMakeFiles/common.dir/ext/glew/glew.c.o"

# External object files for target common
common_EXTERNAL_OBJECTS =

/Users/dnalwqer/Desktop/CG/Dart-Graphics-Assign01-master/bin/mk/libcommon.a: src/common/CMakeFiles/common.dir/image.cpp.o
/Users/dnalwqer/Desktop/CG/Dart-Graphics-Assign01-master/bin/mk/libcommon.a: src/common/CMakeFiles/common.dir/json.cpp.o
/Users/dnalwqer/Desktop/CG/Dart-Graphics-Assign01-master/bin/mk/libcommon.a: src/common/CMakeFiles/common.dir/scene.cpp.o
/Users/dnalwqer/Desktop/CG/Dart-Graphics-Assign01-master/bin/mk/libcommon.a: src/common/CMakeFiles/common.dir/ext/lodepng/lodepng.cpp.o
/Users/dnalwqer/Desktop/CG/Dart-Graphics-Assign01-master/bin/mk/libcommon.a: src/common/CMakeFiles/common.dir/ext/glew/glew.c.o
/Users/dnalwqer/Desktop/CG/Dart-Graphics-Assign01-master/bin/mk/libcommon.a: src/common/CMakeFiles/common.dir/build.make
/Users/dnalwqer/Desktop/CG/Dart-Graphics-Assign01-master/bin/mk/libcommon.a: src/common/CMakeFiles/common.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/dnalwqer/Desktop/CG/Dart-Graphics-Assign01-master-build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Linking CXX static library /Users/dnalwqer/Desktop/CG/Dart-Graphics-Assign01-master/bin/mk/libcommon.a"
	cd /Users/dnalwqer/Desktop/CG/Dart-Graphics-Assign01-master-build/src/common && $(CMAKE_COMMAND) -P CMakeFiles/common.dir/cmake_clean_target.cmake
	cd /Users/dnalwqer/Desktop/CG/Dart-Graphics-Assign01-master-build/src/common && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/common.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
src/common/CMakeFiles/common.dir/build: /Users/dnalwqer/Desktop/CG/Dart-Graphics-Assign01-master/bin/mk/libcommon.a

.PHONY : src/common/CMakeFiles/common.dir/build

src/common/CMakeFiles/common.dir/requires: src/common/CMakeFiles/common.dir/image.cpp.o.requires
src/common/CMakeFiles/common.dir/requires: src/common/CMakeFiles/common.dir/json.cpp.o.requires
src/common/CMakeFiles/common.dir/requires: src/common/CMakeFiles/common.dir/scene.cpp.o.requires
src/common/CMakeFiles/common.dir/requires: src/common/CMakeFiles/common.dir/ext/lodepng/lodepng.cpp.o.requires
src/common/CMakeFiles/common.dir/requires: src/common/CMakeFiles/common.dir/ext/glew/glew.c.o.requires

.PHONY : src/common/CMakeFiles/common.dir/requires

src/common/CMakeFiles/common.dir/clean:
	cd /Users/dnalwqer/Desktop/CG/Dart-Graphics-Assign01-master-build/src/common && $(CMAKE_COMMAND) -P CMakeFiles/common.dir/cmake_clean.cmake
.PHONY : src/common/CMakeFiles/common.dir/clean

src/common/CMakeFiles/common.dir/depend:
	cd /Users/dnalwqer/Desktop/CG/Dart-Graphics-Assign01-master-build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/dnalwqer/Desktop/CG/Dart-Graphics-Assign01-master /Users/dnalwqer/Desktop/CG/Dart-Graphics-Assign01-master/src/common /Users/dnalwqer/Desktop/CG/Dart-Graphics-Assign01-master-build /Users/dnalwqer/Desktop/CG/Dart-Graphics-Assign01-master-build/src/common /Users/dnalwqer/Desktop/CG/Dart-Graphics-Assign01-master-build/src/common/CMakeFiles/common.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/common/CMakeFiles/common.dir/depend
