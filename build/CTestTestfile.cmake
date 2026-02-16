# CMake generated Testfile for 
# Source directory: /Users/manavmadanrawal/Dev/schrodingerssandbox
# Build directory: /Users/manavmadanrawal/Dev/schrodingerssandbox/build
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
include("/Users/manavmadanrawal/Dev/schrodingerssandbox/build/test_constants[1]_include.cmake")
add_test([=[test_special_functions]=] "/Users/manavmadanrawal/Dev/schrodingerssandbox/build/test_special_functions")
set_tests_properties([=[test_special_functions]=] PROPERTIES  _BACKTRACE_TRIPLES "/Users/manavmadanrawal/Dev/schrodingerssandbox/CMakeLists.txt;98;add_test;/Users/manavmadanrawal/Dev/schrodingerssandbox/CMakeLists.txt;0;")
add_test([=[test_hydrogen]=] "/Users/manavmadanrawal/Dev/schrodingerssandbox/build/test_hydrogen")
set_tests_properties([=[test_hydrogen]=] PROPERTIES  _BACKTRACE_TRIPLES "/Users/manavmadanrawal/Dev/schrodingerssandbox/CMakeLists.txt;99;add_test;/Users/manavmadanrawal/Dev/schrodingerssandbox/CMakeLists.txt;0;")
add_test([=[test_slater]=] "/Users/manavmadanrawal/Dev/schrodingerssandbox/build/test_slater")
set_tests_properties([=[test_slater]=] PROPERTIES  _BACKTRACE_TRIPLES "/Users/manavmadanrawal/Dev/schrodingerssandbox/CMakeLists.txt;100;add_test;/Users/manavmadanrawal/Dev/schrodingerssandbox/CMakeLists.txt;0;")
subdirs("external/glfw")
subdirs("external/googletest")
