# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.22

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
CMAKE_SOURCE_DIR = /home/wfz/Yux_Project/TFHE/CBSmode

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/wfz/Yux_Project/TFHE/CBSmode/build

# Include any dependencies generated for this target.
include Yux/CMakeFiles/homoYux2.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include Yux/CMakeFiles/homoYux2.dir/compiler_depend.make

# Include the progress variables for this target.
include Yux/CMakeFiles/homoYux2.dir/progress.make

# Include the compile flags for this target's objects.
include Yux/CMakeFiles/homoYux2.dir/flags.make

Yux/CMakeFiles/homoYux2.dir/Yux2.cpp.o: Yux/CMakeFiles/homoYux2.dir/flags.make
Yux/CMakeFiles/homoYux2.dir/Yux2.cpp.o: ../Yux/Yux2.cpp
Yux/CMakeFiles/homoYux2.dir/Yux2.cpp.o: Yux/CMakeFiles/homoYux2.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/wfz/Yux_Project/TFHE/CBSmode/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object Yux/CMakeFiles/homoYux2.dir/Yux2.cpp.o"
	cd /home/wfz/Yux_Project/TFHE/CBSmode/build/Yux && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT Yux/CMakeFiles/homoYux2.dir/Yux2.cpp.o -MF CMakeFiles/homoYux2.dir/Yux2.cpp.o.d -o CMakeFiles/homoYux2.dir/Yux2.cpp.o -c /home/wfz/Yux_Project/TFHE/CBSmode/Yux/Yux2.cpp

Yux/CMakeFiles/homoYux2.dir/Yux2.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/homoYux2.dir/Yux2.cpp.i"
	cd /home/wfz/Yux_Project/TFHE/CBSmode/build/Yux && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/wfz/Yux_Project/TFHE/CBSmode/Yux/Yux2.cpp > CMakeFiles/homoYux2.dir/Yux2.cpp.i

Yux/CMakeFiles/homoYux2.dir/Yux2.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/homoYux2.dir/Yux2.cpp.s"
	cd /home/wfz/Yux_Project/TFHE/CBSmode/build/Yux && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/wfz/Yux_Project/TFHE/CBSmode/Yux/Yux2.cpp -o CMakeFiles/homoYux2.dir/Yux2.cpp.s

Yux/CMakeFiles/homoYux2.dir/Yux2_homo.cpp.o: Yux/CMakeFiles/homoYux2.dir/flags.make
Yux/CMakeFiles/homoYux2.dir/Yux2_homo.cpp.o: ../Yux/Yux2_homo.cpp
Yux/CMakeFiles/homoYux2.dir/Yux2_homo.cpp.o: Yux/CMakeFiles/homoYux2.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/wfz/Yux_Project/TFHE/CBSmode/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object Yux/CMakeFiles/homoYux2.dir/Yux2_homo.cpp.o"
	cd /home/wfz/Yux_Project/TFHE/CBSmode/build/Yux && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT Yux/CMakeFiles/homoYux2.dir/Yux2_homo.cpp.o -MF CMakeFiles/homoYux2.dir/Yux2_homo.cpp.o.d -o CMakeFiles/homoYux2.dir/Yux2_homo.cpp.o -c /home/wfz/Yux_Project/TFHE/CBSmode/Yux/Yux2_homo.cpp

Yux/CMakeFiles/homoYux2.dir/Yux2_homo.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/homoYux2.dir/Yux2_homo.cpp.i"
	cd /home/wfz/Yux_Project/TFHE/CBSmode/build/Yux && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/wfz/Yux_Project/TFHE/CBSmode/Yux/Yux2_homo.cpp > CMakeFiles/homoYux2.dir/Yux2_homo.cpp.i

Yux/CMakeFiles/homoYux2.dir/Yux2_homo.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/homoYux2.dir/Yux2_homo.cpp.s"
	cd /home/wfz/Yux_Project/TFHE/CBSmode/build/Yux && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/wfz/Yux_Project/TFHE/CBSmode/Yux/Yux2_homo.cpp -o CMakeFiles/homoYux2.dir/Yux2_homo.cpp.s

# Object files for target homoYux2
homoYux2_OBJECTS = \
"CMakeFiles/homoYux2.dir/Yux2.cpp.o" \
"CMakeFiles/homoYux2.dir/Yux2_homo.cpp.o"

# External object files for target homoYux2
homoYux2_EXTERNAL_OBJECTS =

Yux/homoYux2: Yux/CMakeFiles/homoYux2.dir/Yux2.cpp.o
Yux/homoYux2: Yux/CMakeFiles/homoYux2.dir/Yux2_homo.cpp.o
Yux/homoYux2: Yux/CMakeFiles/homoYux2.dir/build.make
Yux/homoYux2: src/libtfhe++.a
Yux/homoYux2: thirdparties/randen/libranden.a
Yux/homoYux2: thirdparties/spqlios/libspqlios.a
Yux/homoYux2: /usr/lib/gcc/x86_64-linux-gnu/11/libgomp.so
Yux/homoYux2: /usr/lib/x86_64-linux-gnu/libpthread.a
Yux/homoYux2: Yux/CMakeFiles/homoYux2.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/wfz/Yux_Project/TFHE/CBSmode/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Linking CXX executable homoYux2"
	cd /home/wfz/Yux_Project/TFHE/CBSmode/build/Yux && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/homoYux2.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
Yux/CMakeFiles/homoYux2.dir/build: Yux/homoYux2
.PHONY : Yux/CMakeFiles/homoYux2.dir/build

Yux/CMakeFiles/homoYux2.dir/clean:
	cd /home/wfz/Yux_Project/TFHE/CBSmode/build/Yux && $(CMAKE_COMMAND) -P CMakeFiles/homoYux2.dir/cmake_clean.cmake
.PHONY : Yux/CMakeFiles/homoYux2.dir/clean

Yux/CMakeFiles/homoYux2.dir/depend:
	cd /home/wfz/Yux_Project/TFHE/CBSmode/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/wfz/Yux_Project/TFHE/CBSmode /home/wfz/Yux_Project/TFHE/CBSmode/Yux /home/wfz/Yux_Project/TFHE/CBSmode/build /home/wfz/Yux_Project/TFHE/CBSmode/build/Yux /home/wfz/Yux_Project/TFHE/CBSmode/build/Yux/CMakeFiles/homoYux2.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : Yux/CMakeFiles/homoYux2.dir/depend

