# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.28

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
CMAKE_COMMAND = /opt/homebrew/Cellar/cmake/3.28.3/bin/cmake

# The command to remove a file.
RM = /opt/homebrew/Cellar/cmake/3.28.3/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/jonathanchen/Desktop/pbrt-v3

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/jonathanchen/Desktop/pbrt-v3/build

# Include any dependencies generated for this target.
include src/ext/ptex/src/utils/CMakeFiles/ptxinfo.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include src/ext/ptex/src/utils/CMakeFiles/ptxinfo.dir/compiler_depend.make

# Include the progress variables for this target.
include src/ext/ptex/src/utils/CMakeFiles/ptxinfo.dir/progress.make

# Include the compile flags for this target's objects.
include src/ext/ptex/src/utils/CMakeFiles/ptxinfo.dir/flags.make

src/ext/ptex/src/utils/CMakeFiles/ptxinfo.dir/ptxinfo.cpp.o: src/ext/ptex/src/utils/CMakeFiles/ptxinfo.dir/flags.make
src/ext/ptex/src/utils/CMakeFiles/ptxinfo.dir/ptxinfo.cpp.o: /Users/jonathanchen/Desktop/pbrt-v3/src/ext/ptex/src/utils/ptxinfo.cpp
src/ext/ptex/src/utils/CMakeFiles/ptxinfo.dir/ptxinfo.cpp.o: src/ext/ptex/src/utils/CMakeFiles/ptxinfo.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/Users/jonathanchen/Desktop/pbrt-v3/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object src/ext/ptex/src/utils/CMakeFiles/ptxinfo.dir/ptxinfo.cpp.o"
	cd /Users/jonathanchen/Desktop/pbrt-v3/build/src/ext/ptex/src/utils && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT src/ext/ptex/src/utils/CMakeFiles/ptxinfo.dir/ptxinfo.cpp.o -MF CMakeFiles/ptxinfo.dir/ptxinfo.cpp.o.d -o CMakeFiles/ptxinfo.dir/ptxinfo.cpp.o -c /Users/jonathanchen/Desktop/pbrt-v3/src/ext/ptex/src/utils/ptxinfo.cpp

src/ext/ptex/src/utils/CMakeFiles/ptxinfo.dir/ptxinfo.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/ptxinfo.dir/ptxinfo.cpp.i"
	cd /Users/jonathanchen/Desktop/pbrt-v3/build/src/ext/ptex/src/utils && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/jonathanchen/Desktop/pbrt-v3/src/ext/ptex/src/utils/ptxinfo.cpp > CMakeFiles/ptxinfo.dir/ptxinfo.cpp.i

src/ext/ptex/src/utils/CMakeFiles/ptxinfo.dir/ptxinfo.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/ptxinfo.dir/ptxinfo.cpp.s"
	cd /Users/jonathanchen/Desktop/pbrt-v3/build/src/ext/ptex/src/utils && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/jonathanchen/Desktop/pbrt-v3/src/ext/ptex/src/utils/ptxinfo.cpp -o CMakeFiles/ptxinfo.dir/ptxinfo.cpp.s

# Object files for target ptxinfo
ptxinfo_OBJECTS = \
"CMakeFiles/ptxinfo.dir/ptxinfo.cpp.o"

# External object files for target ptxinfo
ptxinfo_EXTERNAL_OBJECTS =

src/ext/ptex/src/utils/ptxinfo: src/ext/ptex/src/utils/CMakeFiles/ptxinfo.dir/ptxinfo.cpp.o
src/ext/ptex/src/utils/ptxinfo: src/ext/ptex/src/utils/CMakeFiles/ptxinfo.dir/build.make
src/ext/ptex/src/utils/ptxinfo: src/ext/ptex/src/ptex/libPtex.a
src/ext/ptex/src/utils/ptxinfo: /Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX14.0.sdk/usr/lib/libz.tbd
src/ext/ptex/src/utils/ptxinfo: src/ext/ptex/src/utils/CMakeFiles/ptxinfo.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --bold --progress-dir=/Users/jonathanchen/Desktop/pbrt-v3/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable ptxinfo"
	cd /Users/jonathanchen/Desktop/pbrt-v3/build/src/ext/ptex/src/utils && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/ptxinfo.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
src/ext/ptex/src/utils/CMakeFiles/ptxinfo.dir/build: src/ext/ptex/src/utils/ptxinfo
.PHONY : src/ext/ptex/src/utils/CMakeFiles/ptxinfo.dir/build

src/ext/ptex/src/utils/CMakeFiles/ptxinfo.dir/clean:
	cd /Users/jonathanchen/Desktop/pbrt-v3/build/src/ext/ptex/src/utils && $(CMAKE_COMMAND) -P CMakeFiles/ptxinfo.dir/cmake_clean.cmake
.PHONY : src/ext/ptex/src/utils/CMakeFiles/ptxinfo.dir/clean

src/ext/ptex/src/utils/CMakeFiles/ptxinfo.dir/depend:
	cd /Users/jonathanchen/Desktop/pbrt-v3/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/jonathanchen/Desktop/pbrt-v3 /Users/jonathanchen/Desktop/pbrt-v3/src/ext/ptex/src/utils /Users/jonathanchen/Desktop/pbrt-v3/build /Users/jonathanchen/Desktop/pbrt-v3/build/src/ext/ptex/src/utils /Users/jonathanchen/Desktop/pbrt-v3/build/src/ext/ptex/src/utils/CMakeFiles/ptxinfo.dir/DependInfo.cmake "--color=$(COLOR)"
.PHONY : src/ext/ptex/src/utils/CMakeFiles/ptxinfo.dir/depend

