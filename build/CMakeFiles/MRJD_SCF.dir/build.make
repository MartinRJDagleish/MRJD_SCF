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
include CMakeFiles/MRJD_SCF.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/MRJD_SCF.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/MRJD_SCF.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/MRJD_SCF.dir/flags.make

CMakeFiles/MRJD_SCF.dir/main.cpp.o: CMakeFiles/MRJD_SCF.dir/flags.make
CMakeFiles/MRJD_SCF.dir/main.cpp.o: /run/media/mrjd/DATA/Programming/C++/Crawford_Grp_MRJD/main.cpp
CMakeFiles/MRJD_SCF.dir/main.cpp.o: CMakeFiles/MRJD_SCF.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/run/media/mrjd/DATA/Programming/C++/Crawford_Grp_MRJD/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/MRJD_SCF.dir/main.cpp.o"
	/usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/MRJD_SCF.dir/main.cpp.o -MF CMakeFiles/MRJD_SCF.dir/main.cpp.o.d -o CMakeFiles/MRJD_SCF.dir/main.cpp.o -c /run/media/mrjd/DATA/Programming/C++/Crawford_Grp_MRJD/main.cpp

CMakeFiles/MRJD_SCF.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/MRJD_SCF.dir/main.cpp.i"
	/usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /run/media/mrjd/DATA/Programming/C++/Crawford_Grp_MRJD/main.cpp > CMakeFiles/MRJD_SCF.dir/main.cpp.i

CMakeFiles/MRJD_SCF.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/MRJD_SCF.dir/main.cpp.s"
	/usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /run/media/mrjd/DATA/Programming/C++/Crawford_Grp_MRJD/main.cpp -o CMakeFiles/MRJD_SCF.dir/main.cpp.s

# Object files for target MRJD_SCF
MRJD_SCF_OBJECTS = \
"CMakeFiles/MRJD_SCF.dir/main.cpp.o"

# External object files for target MRJD_SCF
MRJD_SCF_EXTERNAL_OBJECTS =

MRJD_SCF: CMakeFiles/MRJD_SCF.dir/main.cpp.o
MRJD_SCF: CMakeFiles/MRJD_SCF.dir/build.make
MRJD_SCF: CMakeFiles/MRJD_SCF.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/run/media/mrjd/DATA/Programming/C++/Crawford_Grp_MRJD/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable MRJD_SCF"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/MRJD_SCF.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/MRJD_SCF.dir/build: MRJD_SCF
.PHONY : CMakeFiles/MRJD_SCF.dir/build

CMakeFiles/MRJD_SCF.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/MRJD_SCF.dir/cmake_clean.cmake
.PHONY : CMakeFiles/MRJD_SCF.dir/clean

CMakeFiles/MRJD_SCF.dir/depend:
	cd /run/media/mrjd/DATA/Programming/C++/Crawford_Grp_MRJD/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /run/media/mrjd/DATA/Programming/C++/Crawford_Grp_MRJD /run/media/mrjd/DATA/Programming/C++/Crawford_Grp_MRJD /run/media/mrjd/DATA/Programming/C++/Crawford_Grp_MRJD/build /run/media/mrjd/DATA/Programming/C++/Crawford_Grp_MRJD/build /run/media/mrjd/DATA/Programming/C++/Crawford_Grp_MRJD/build/CMakeFiles/MRJD_SCF.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/MRJD_SCF.dir/depend

