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
include src/ext/openexr/IlmBase/Iex/CMakeFiles/Iex.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include src/ext/openexr/IlmBase/Iex/CMakeFiles/Iex.dir/compiler_depend.make

# Include the progress variables for this target.
include src/ext/openexr/IlmBase/Iex/CMakeFiles/Iex.dir/progress.make

# Include the compile flags for this target's objects.
include src/ext/openexr/IlmBase/Iex/CMakeFiles/Iex.dir/flags.make

src/ext/openexr/IlmBase/Iex/CMakeFiles/Iex.dir/IexBaseExc.cpp.o: src/ext/openexr/IlmBase/Iex/CMakeFiles/Iex.dir/flags.make
src/ext/openexr/IlmBase/Iex/CMakeFiles/Iex.dir/IexBaseExc.cpp.o: /Users/jonathanchen/Desktop/pbrt-v3/src/ext/openexr/IlmBase/Iex/IexBaseExc.cpp
src/ext/openexr/IlmBase/Iex/CMakeFiles/Iex.dir/IexBaseExc.cpp.o: src/ext/openexr/IlmBase/Iex/CMakeFiles/Iex.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/Users/jonathanchen/Desktop/pbrt-v3/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object src/ext/openexr/IlmBase/Iex/CMakeFiles/Iex.dir/IexBaseExc.cpp.o"
	cd /Users/jonathanchen/Desktop/pbrt-v3/build/src/ext/openexr/IlmBase/Iex && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT src/ext/openexr/IlmBase/Iex/CMakeFiles/Iex.dir/IexBaseExc.cpp.o -MF CMakeFiles/Iex.dir/IexBaseExc.cpp.o.d -o CMakeFiles/Iex.dir/IexBaseExc.cpp.o -c /Users/jonathanchen/Desktop/pbrt-v3/src/ext/openexr/IlmBase/Iex/IexBaseExc.cpp

src/ext/openexr/IlmBase/Iex/CMakeFiles/Iex.dir/IexBaseExc.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/Iex.dir/IexBaseExc.cpp.i"
	cd /Users/jonathanchen/Desktop/pbrt-v3/build/src/ext/openexr/IlmBase/Iex && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/jonathanchen/Desktop/pbrt-v3/src/ext/openexr/IlmBase/Iex/IexBaseExc.cpp > CMakeFiles/Iex.dir/IexBaseExc.cpp.i

src/ext/openexr/IlmBase/Iex/CMakeFiles/Iex.dir/IexBaseExc.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/Iex.dir/IexBaseExc.cpp.s"
	cd /Users/jonathanchen/Desktop/pbrt-v3/build/src/ext/openexr/IlmBase/Iex && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/jonathanchen/Desktop/pbrt-v3/src/ext/openexr/IlmBase/Iex/IexBaseExc.cpp -o CMakeFiles/Iex.dir/IexBaseExc.cpp.s

src/ext/openexr/IlmBase/Iex/CMakeFiles/Iex.dir/IexThrowErrnoExc.cpp.o: src/ext/openexr/IlmBase/Iex/CMakeFiles/Iex.dir/flags.make
src/ext/openexr/IlmBase/Iex/CMakeFiles/Iex.dir/IexThrowErrnoExc.cpp.o: /Users/jonathanchen/Desktop/pbrt-v3/src/ext/openexr/IlmBase/Iex/IexThrowErrnoExc.cpp
src/ext/openexr/IlmBase/Iex/CMakeFiles/Iex.dir/IexThrowErrnoExc.cpp.o: src/ext/openexr/IlmBase/Iex/CMakeFiles/Iex.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/Users/jonathanchen/Desktop/pbrt-v3/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object src/ext/openexr/IlmBase/Iex/CMakeFiles/Iex.dir/IexThrowErrnoExc.cpp.o"
	cd /Users/jonathanchen/Desktop/pbrt-v3/build/src/ext/openexr/IlmBase/Iex && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT src/ext/openexr/IlmBase/Iex/CMakeFiles/Iex.dir/IexThrowErrnoExc.cpp.o -MF CMakeFiles/Iex.dir/IexThrowErrnoExc.cpp.o.d -o CMakeFiles/Iex.dir/IexThrowErrnoExc.cpp.o -c /Users/jonathanchen/Desktop/pbrt-v3/src/ext/openexr/IlmBase/Iex/IexThrowErrnoExc.cpp

src/ext/openexr/IlmBase/Iex/CMakeFiles/Iex.dir/IexThrowErrnoExc.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/Iex.dir/IexThrowErrnoExc.cpp.i"
	cd /Users/jonathanchen/Desktop/pbrt-v3/build/src/ext/openexr/IlmBase/Iex && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/jonathanchen/Desktop/pbrt-v3/src/ext/openexr/IlmBase/Iex/IexThrowErrnoExc.cpp > CMakeFiles/Iex.dir/IexThrowErrnoExc.cpp.i

src/ext/openexr/IlmBase/Iex/CMakeFiles/Iex.dir/IexThrowErrnoExc.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/Iex.dir/IexThrowErrnoExc.cpp.s"
	cd /Users/jonathanchen/Desktop/pbrt-v3/build/src/ext/openexr/IlmBase/Iex && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/jonathanchen/Desktop/pbrt-v3/src/ext/openexr/IlmBase/Iex/IexThrowErrnoExc.cpp -o CMakeFiles/Iex.dir/IexThrowErrnoExc.cpp.s

# Object files for target Iex
Iex_OBJECTS = \
"CMakeFiles/Iex.dir/IexBaseExc.cpp.o" \
"CMakeFiles/Iex.dir/IexThrowErrnoExc.cpp.o"

# External object files for target Iex
Iex_EXTERNAL_OBJECTS =

src/ext/openexr/IlmBase/Iex/libIex.a: src/ext/openexr/IlmBase/Iex/CMakeFiles/Iex.dir/IexBaseExc.cpp.o
src/ext/openexr/IlmBase/Iex/libIex.a: src/ext/openexr/IlmBase/Iex/CMakeFiles/Iex.dir/IexThrowErrnoExc.cpp.o
src/ext/openexr/IlmBase/Iex/libIex.a: src/ext/openexr/IlmBase/Iex/CMakeFiles/Iex.dir/build.make
src/ext/openexr/IlmBase/Iex/libIex.a: src/ext/openexr/IlmBase/Iex/CMakeFiles/Iex.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --bold --progress-dir=/Users/jonathanchen/Desktop/pbrt-v3/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Linking CXX static library libIex.a"
	cd /Users/jonathanchen/Desktop/pbrt-v3/build/src/ext/openexr/IlmBase/Iex && $(CMAKE_COMMAND) -P CMakeFiles/Iex.dir/cmake_clean_target.cmake
	cd /Users/jonathanchen/Desktop/pbrt-v3/build/src/ext/openexr/IlmBase/Iex && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/Iex.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
src/ext/openexr/IlmBase/Iex/CMakeFiles/Iex.dir/build: src/ext/openexr/IlmBase/Iex/libIex.a
.PHONY : src/ext/openexr/IlmBase/Iex/CMakeFiles/Iex.dir/build

src/ext/openexr/IlmBase/Iex/CMakeFiles/Iex.dir/clean:
	cd /Users/jonathanchen/Desktop/pbrt-v3/build/src/ext/openexr/IlmBase/Iex && $(CMAKE_COMMAND) -P CMakeFiles/Iex.dir/cmake_clean.cmake
.PHONY : src/ext/openexr/IlmBase/Iex/CMakeFiles/Iex.dir/clean

src/ext/openexr/IlmBase/Iex/CMakeFiles/Iex.dir/depend:
	cd /Users/jonathanchen/Desktop/pbrt-v3/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/jonathanchen/Desktop/pbrt-v3 /Users/jonathanchen/Desktop/pbrt-v3/src/ext/openexr/IlmBase/Iex /Users/jonathanchen/Desktop/pbrt-v3/build /Users/jonathanchen/Desktop/pbrt-v3/build/src/ext/openexr/IlmBase/Iex /Users/jonathanchen/Desktop/pbrt-v3/build/src/ext/openexr/IlmBase/Iex/CMakeFiles/Iex.dir/DependInfo.cmake "--color=$(COLOR)"
.PHONY : src/ext/openexr/IlmBase/Iex/CMakeFiles/Iex.dir/depend

