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
CMAKE_COMMAND = /opt/clion-2018.1.1/bin/cmake/bin/cmake

# The command to remove a file.
RM = /opt/clion-2018.1.1/bin/cmake/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/wt/Project/bwa-HashMEM

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/wt/Project/bwa-HashMEM/cmake-build-debug

# Include any dependencies generated for this target.
include CMakeFiles/bwa_HashMEM.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/bwa_HashMEM.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/bwa_HashMEM.dir/flags.make

CMakeFiles/bwa_HashMEM.dir/main.c.o: CMakeFiles/bwa_HashMEM.dir/flags.make
CMakeFiles/bwa_HashMEM.dir/main.c.o: ../main.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/wt/Project/bwa-HashMEM/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building C object CMakeFiles/bwa_HashMEM.dir/main.c.o"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/bwa_HashMEM.dir/main.c.o   -c /home/wt/Project/bwa-HashMEM/main.c

CMakeFiles/bwa_HashMEM.dir/main.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/bwa_HashMEM.dir/main.c.i"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/wt/Project/bwa-HashMEM/main.c > CMakeFiles/bwa_HashMEM.dir/main.c.i

CMakeFiles/bwa_HashMEM.dir/main.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/bwa_HashMEM.dir/main.c.s"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/wt/Project/bwa-HashMEM/main.c -o CMakeFiles/bwa_HashMEM.dir/main.c.s

CMakeFiles/bwa_HashMEM.dir/main.c.o.requires:

.PHONY : CMakeFiles/bwa_HashMEM.dir/main.c.o.requires

CMakeFiles/bwa_HashMEM.dir/main.c.o.provides: CMakeFiles/bwa_HashMEM.dir/main.c.o.requires
	$(MAKE) -f CMakeFiles/bwa_HashMEM.dir/build.make CMakeFiles/bwa_HashMEM.dir/main.c.o.provides.build
.PHONY : CMakeFiles/bwa_HashMEM.dir/main.c.o.provides

CMakeFiles/bwa_HashMEM.dir/main.c.o.provides.build: CMakeFiles/bwa_HashMEM.dir/main.c.o


# Object files for target bwa_HashMEM
bwa_HashMEM_OBJECTS = \
"CMakeFiles/bwa_HashMEM.dir/main.c.o"

# External object files for target bwa_HashMEM
bwa_HashMEM_EXTERNAL_OBJECTS =

bwa_HashMEM: CMakeFiles/bwa_HashMEM.dir/main.c.o
bwa_HashMEM: CMakeFiles/bwa_HashMEM.dir/build.make
bwa_HashMEM: CMakeFiles/bwa_HashMEM.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/wt/Project/bwa-HashMEM/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking C executable bwa_HashMEM"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/bwa_HashMEM.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/bwa_HashMEM.dir/build: bwa_HashMEM

.PHONY : CMakeFiles/bwa_HashMEM.dir/build

CMakeFiles/bwa_HashMEM.dir/requires: CMakeFiles/bwa_HashMEM.dir/main.c.o.requires

.PHONY : CMakeFiles/bwa_HashMEM.dir/requires

CMakeFiles/bwa_HashMEM.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/bwa_HashMEM.dir/cmake_clean.cmake
.PHONY : CMakeFiles/bwa_HashMEM.dir/clean

CMakeFiles/bwa_HashMEM.dir/depend:
	cd /home/wt/Project/bwa-HashMEM/cmake-build-debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/wt/Project/bwa-HashMEM /home/wt/Project/bwa-HashMEM /home/wt/Project/bwa-HashMEM/cmake-build-debug /home/wt/Project/bwa-HashMEM/cmake-build-debug /home/wt/Project/bwa-HashMEM/cmake-build-debug/CMakeFiles/bwa_HashMEM.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/bwa_HashMEM.dir/depend

