# CMake generated Testfile for 
# Source directory: /Users/jonathanchen/Desktop/pbrt-v3/src/ext/ptex/src/tests
# Build directory: /Users/jonathanchen/Desktop/pbrt-v3/build/src/ext/ptex/src/tests
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(wtest "/Users/jonathanchen/Desktop/pbrt-v3/build/src/ext/ptex/src/tests/wtest")
set_tests_properties(wtest PROPERTIES  _BACKTRACE_TRIPLES "/Users/jonathanchen/Desktop/pbrt-v3/src/ext/ptex/src/tests/CMakeLists.txt;32;add_test;/Users/jonathanchen/Desktop/pbrt-v3/src/ext/ptex/src/tests/CMakeLists.txt;0;")
add_test(rtest "/opt/homebrew/Cellar/cmake/3.28.3/bin/cmake" "-DOUT=/Users/jonathanchen/Desktop/pbrt-v3/build/src/ext/ptex/src/tests/rtest.out" "-DDATA=/Users/jonathanchen/Desktop/pbrt-v3/src/ext/ptex/src/tests/rtestok.dat" "-DCMD=/Users/jonathanchen/Desktop/pbrt-v3/build/src/ext/ptex/src/tests/rtest" "-P" "/Users/jonathanchen/Desktop/pbrt-v3/src/ext/ptex/src/tests/compare_test.cmake")
set_tests_properties(rtest PROPERTIES  _BACKTRACE_TRIPLES "/Users/jonathanchen/Desktop/pbrt-v3/src/ext/ptex/src/tests/CMakeLists.txt;23;add_test;/Users/jonathanchen/Desktop/pbrt-v3/src/ext/ptex/src/tests/CMakeLists.txt;33;add_compare_test;/Users/jonathanchen/Desktop/pbrt-v3/src/ext/ptex/src/tests/CMakeLists.txt;0;")
add_test(ftest "/opt/homebrew/Cellar/cmake/3.28.3/bin/cmake" "-DOUT=/Users/jonathanchen/Desktop/pbrt-v3/build/src/ext/ptex/src/tests/ftest.out" "-DDATA=/Users/jonathanchen/Desktop/pbrt-v3/src/ext/ptex/src/tests/ftestok.dat" "-DCMD=/Users/jonathanchen/Desktop/pbrt-v3/build/src/ext/ptex/src/tests/ftest" "-P" "/Users/jonathanchen/Desktop/pbrt-v3/src/ext/ptex/src/tests/compare_test.cmake")
set_tests_properties(ftest PROPERTIES  _BACKTRACE_TRIPLES "/Users/jonathanchen/Desktop/pbrt-v3/src/ext/ptex/src/tests/CMakeLists.txt;23;add_test;/Users/jonathanchen/Desktop/pbrt-v3/src/ext/ptex/src/tests/CMakeLists.txt;34;add_compare_test;/Users/jonathanchen/Desktop/pbrt-v3/src/ext/ptex/src/tests/CMakeLists.txt;0;")
add_test(halftest "/Users/jonathanchen/Desktop/pbrt-v3/build/src/ext/ptex/src/tests/halftest")
set_tests_properties(halftest PROPERTIES  _BACKTRACE_TRIPLES "/Users/jonathanchen/Desktop/pbrt-v3/src/ext/ptex/src/tests/CMakeLists.txt;35;add_test;/Users/jonathanchen/Desktop/pbrt-v3/src/ext/ptex/src/tests/CMakeLists.txt;0;")
