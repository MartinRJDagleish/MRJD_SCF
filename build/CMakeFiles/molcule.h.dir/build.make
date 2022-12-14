# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.24

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /run/media/mrjd/DATA/Programming/C++/Crawford_Grp_MRJD

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /run/media/mrjd/DATA/Programming/C++/Crawford_Grp_MRJD/build

# Include any dependencies generated for this target.
include CMakeFiles/molcule.h.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/molcule.h.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/molcule.h.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/molcule.h.dir/flags.make

CMakeFiles/molcule.h.dir/molecule.cc.o: CMakeFiles/molcule.h.dir/flags.make
CMakeFiles/molcule.h.dir/molecule.cc.o: /run/media/mrjd/DATA/Programming/C++/Crawford_Grp_MRJD/molecule.cc
CMakeFiles/molcule.h.dir/molecule.cc.o: CMakeFiles/molcule.h.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/run/media/mrjd/DATA/Programming/C++/Crawford_Grp_MRJD/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/molcule.h.dir/molecule.cc.o"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/molcule.h.dir/molecule.cc.o -MF CMakeFiles/molcule.h.dir/molecule.cc.o.d -o CMakeFiles/molcule.h.dir/molecule.cc.o -c /run/media/mrjd/DATA/Programming/C++/Crawford_Grp_MRJD/molecule.cc

CMakeFiles/molcule.h.dir/molecule.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/molcule.h.dir/molecule.cc.i"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /run/media/mrjd/DATA/Programming/C++/Crawford_Grp_MRJD/molecule.cc > CMakeFiles/molcule.h.dir/molecule.cc.i

CMakeFiles/molcule.h.dir/molecule.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/molcule.h.dir/molecule.cc.s"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /run/media/mrjd/DATA/Programming/C++/Crawford_Grp_MRJD/molecule.cc -o CMakeFiles/molcule.h.dir/molecule.cc.s

# Object files for target molcule.h
molcule_h_OBJECTS = \
"CMakeFiles/molcule.h.dir/molecule.cc.o"

# External object files for target molcule.h
molcule_h_EXTERNAL_OBJECTS =

libmolcule.h.a: CMakeFiles/molcule.h.dir/molecule.cc.o
libmolcule.h.a: CMakeFiles/molcule.h.dir/build.make
libmolcule.h.a: CMakeFiles/molcule.h.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/run/media/mrjd/DATA/Programming/C++/Crawford_Grp_MRJD/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX static library libmolcule.h.a"
	$(CMAKE_COMMAND) -P CMakeFiles/molcule.h.dir/cmake_clean_target.cmake
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/molcule.h.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/molcule.h.dir/build: libmolcule.h.a
.PHONY : CMakeFiles/molcule.h.dir/build

CMakeFiles/molcule.h.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/molcule.h.dir/cmake_clean.cmake
.PHONY : CMakeFiles/molcule.h.dir/clean

CMakeFiles/molcule.h.dir/depend:
	cd /run/media/mrjd/DATA/Programming/C++/Crawford_Grp_MRJD/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /run/media/mrjd/DATA/Programming/C++/Crawford_Grp_MRJD /run/media/mrjd/DATA/Programming/C++/Crawford_Grp_MRJD /run/media/mrjd/DATA/Programming/C++/Crawford_Grp_MRJD/build /run/media/mrjd/DATA/Programming/C++/Crawford_Grp_MRJD/build /run/media/mrjd/DATA/Programming/C++/Crawford_Grp_MRJD/build/CMakeFiles/molcule.h.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/molcule.h.dir/depend

