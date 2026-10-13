cmake_minimum_required(VERSION 3.16)
project(marsq_halfmesh_regression NONE)
enable_testing()
find_package(Python3 REQUIRED COMPONENTS Interpreter)
add_test(NAME halfmesh_uninitialised_profiles
  COMMAND ${Python3_EXECUTABLE} -m unittest discover -s "${CMAKE_CURRENT_SOURCE_DIR}" -p "test_*.py" -v)
set_tests_properties(halfmesh_uninitialised_profiles PROPERTIES TIMEOUT 900
  ENVIRONMENT "TMPDIR=${CMAKE_CURRENT_BINARY_DIR}")
