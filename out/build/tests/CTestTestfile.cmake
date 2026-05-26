# CMake generated Testfile for 
# Source directory: D:/GitHub/1DNavierStokes-Transient/tests
# Build directory: D:/GitHub/1DNavierStokes-Transient/out/build/tests
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test([=[tsensor_test]=] "tsensor" "--test")
set_tests_properties([=[tsensor_test]=] PROPERTIES  _BACKTRACE_TRIPLES "D:/GitHub/1DNavierStokes-Transient/tests/CMakeLists.txt;20;add_test;D:/GitHub/1DNavierStokes-Transient/tests/CMakeLists.txt;0;")
subdirs("../_deps/googletest-build")
