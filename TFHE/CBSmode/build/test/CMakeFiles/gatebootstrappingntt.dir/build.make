# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.30

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
CMAKE_COMMAND = /usr/local/bin/cmake

# The command to remove a file.
RM = /usr/local/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/wangfangzhen/Yux/CBSmode

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/wangfangzhen/Yux/CBSmode/build

# Include any dependencies generated for this target.
include test/CMakeFiles/gatebootstrappingntt.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include test/CMakeFiles/gatebootstrappingntt.dir/compiler_depend.make

# Include the progress variables for this target.
include test/CMakeFiles/gatebootstrappingntt.dir/progress.make

# Include the compile flags for this target's objects.
include test/CMakeFiles/gatebootstrappingntt.dir/flags.make

test/CMakeFiles/gatebootstrappingntt.dir/gatebootstrappingntt.cpp.o: test/CMakeFiles/gatebootstrappingntt.dir/flags.make
test/CMakeFiles/gatebootstrappingntt.dir/gatebootstrappingntt.cpp.o: /home/wangfangzhen/Yux/CBSmode/test/gatebootstrappingntt.cpp
test/CMakeFiles/gatebootstrappingntt.dir/gatebootstrappingntt.cpp.o: test/CMakeFiles/gatebootstrappingntt.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/wangfangzhen/Yux/CBSmode/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object test/CMakeFiles/gatebootstrappingntt.dir/gatebootstrappingntt.cpp.o"
	cd /home/wangfangzhen/Yux/CBSmode/build/test && /usr/local/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT test/CMakeFiles/gatebootstrappingntt.dir/gatebootstrappingntt.cpp.o -MF CMakeFiles/gatebootstrappingntt.dir/gatebootstrappingntt.cpp.o.d -o CMakeFiles/gatebootstrappingntt.dir/gatebootstrappingntt.cpp.o -c /home/wangfangzhen/Yux/CBSmode/test/gatebootstrappingntt.cpp

test/CMakeFiles/gatebootstrappingntt.dir/gatebootstrappingntt.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/gatebootstrappingntt.dir/gatebootstrappingntt.cpp.i"
	cd /home/wangfangzhen/Yux/CBSmode/build/test && /usr/local/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/wangfangzhen/Yux/CBSmode/test/gatebootstrappingntt.cpp > CMakeFiles/gatebootstrappingntt.dir/gatebootstrappingntt.cpp.i

test/CMakeFiles/gatebootstrappingntt.dir/gatebootstrappingntt.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/gatebootstrappingntt.dir/gatebootstrappingntt.cpp.s"
	cd /home/wangfangzhen/Yux/CBSmode/build/test && /usr/local/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/wangfangzhen/Yux/CBSmode/test/gatebootstrappingntt.cpp -o CMakeFiles/gatebootstrappingntt.dir/gatebootstrappingntt.cpp.s

# Object files for target gatebootstrappingntt
gatebootstrappingntt_OBJECTS = \
"CMakeFiles/gatebootstrappingntt.dir/gatebootstrappingntt.cpp.o"

# External object files for target gatebootstrappingntt
gatebootstrappingntt_EXTERNAL_OBJECTS =

test/gatebootstrappingntt: test/CMakeFiles/gatebootstrappingntt.dir/gatebootstrappingntt.cpp.o
test/gatebootstrappingntt: test/CMakeFiles/gatebootstrappingntt.dir/build.make
test/gatebootstrappingntt: src/libtfhe++.a
test/gatebootstrappingntt: thirdparties/randen/libranden.a
test/gatebootstrappingntt: thirdparties/spqlios/libspqlios.a
test/gatebootstrappingntt: test/CMakeFiles/gatebootstrappingntt.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --bold --progress-dir=/home/wangfangzhen/Yux/CBSmode/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable gatebootstrappingntt"
	cd /home/wangfangzhen/Yux/CBSmode/build/test && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/gatebootstrappingntt.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
test/CMakeFiles/gatebootstrappingntt.dir/build: test/gatebootstrappingntt
.PHONY : test/CMakeFiles/gatebootstrappingntt.dir/build

test/CMakeFiles/gatebootstrappingntt.dir/clean:
	cd /home/wangfangzhen/Yux/CBSmode/build/test && $(CMAKE_COMMAND) -P CMakeFiles/gatebootstrappingntt.dir/cmake_clean.cmake
.PHONY : test/CMakeFiles/gatebootstrappingntt.dir/clean

test/CMakeFiles/gatebootstrappingntt.dir/depend:
	cd /home/wangfangzhen/Yux/CBSmode/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/wangfangzhen/Yux/CBSmode /home/wangfangzhen/Yux/CBSmode/test /home/wangfangzhen/Yux/CBSmode/build /home/wangfangzhen/Yux/CBSmode/build/test /home/wangfangzhen/Yux/CBSmode/build/test/CMakeFiles/gatebootstrappingntt.dir/DependInfo.cmake "--color=$(COLOR)"
.PHONY : test/CMakeFiles/gatebootstrappingntt.dir/depend
