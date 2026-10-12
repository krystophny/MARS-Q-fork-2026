cmake_minimum_required(VERSION 3.16)
project(marsq_chease_regression NONE)
enable_testing()
find_package(Python3 REQUIRED COMPONENTS Interpreter)
add_test(NAME chease_regression
  COMMAND ${Python3_EXECUTABLE} -m unittest discover -s "${CMAKE_CURRENT_SOURCE_DIR}" -p "test_*.py" -v)
set_tests_properties(chease_regression PROPERTIES TIMEOUT 300
  ENVIRONMENT "TMPDIR=${CMAKE_CURRENT_BINARY_DIR}")
