cmake_minimum_required(VERSION 3.20)
project(chease_profile_integral_oracle LANGUAGES Fortran)
enable_testing()
find_package(Python3 REQUIRED COMPONENTS Interpreter)
set(CHEASE_ROOT "${CMAKE_CURRENT_LIST_DIR}/../../CheaseMerge" CACHE PATH "Pinned native source")
execute_process(COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_LIST_DIR}/extract.py"
    "${CHEASE_ROOT}/chease.f" "${CMAKE_CURRENT_BINARY_DIR}" COMMAND_ERROR_IS_FATAL ANY)
add_executable(native_profile_integral "${CMAKE_CURRENT_BINARY_DIR}/native_extracted.f" oracle.f90)
target_include_directories(native_profile_integral PRIVATE "${CHEASE_ROOT}")
target_compile_options(native_profile_integral PRIVATE -fdefault-real-8 -fdefault-double-8
    -ffixed-line-length-none -std=legacy)
foreach(degree RANGE 0 3)
    foreach(direction RANGE 0 1)
        add_test(NAME polynomial_${degree}_direction_${direction}
            COMMAND native_profile_integral ${degree} ${direction})
    endforeach()
endforeach()
foreach(degree RANGE 2 3)
    foreach(direction RANGE 0 1)
        add_test(NAME reject_lost_curvature_${degree}_${direction}
            COMMAND native_profile_integral ${degree} ${direction} mutate)
        set_tests_properties(reject_lost_curvature_${degree}_${direction} PROPERTIES WILL_FAIL TRUE)
    endforeach()
endforeach()
