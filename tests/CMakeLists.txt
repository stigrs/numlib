cmake_minimum_required(VERSION 3.4.3)
project(numlib_tests CXX)

include_directories(${numlib_SOURCE_DIR}/include)
if(NOT APPLE)
	include_directories(${BLAS_INCLUDE_DIR})
	include_directories(${LAPACKE_INCLUDE_DIR})
endif()

# Copy *.inp files to build directory.
file(GLOB INP_FILES RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}/*.inp")

foreach(inp_file_name ${INP_FILES})
    file(COPY "${inp_file_name}" DESTINATION "${CMAKE_BINARY_DIR}/tests")
endforeach()

# Prepare "Catch" library for other executables.
set(CATCH_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR})
add_library(Catch INTERFACE)
target_include_directories(Catch INTERFACE ${CATCH_INCLUDE_DIR})

# Create the main executable for each test.
add_library(test_catch STATIC test.cpp)
target_link_libraries(
    test_catch
    Catch
)

function(add_numlib_test name)
    add_executable(${name} ${name}.cpp)
    target_link_libraries (
        ${program} 
		${BLAS_LIBRARIES}
		test_catch
    ) 
    add_dependencies(${name} Catch)
    add_test(NAME ${name} COMMAND ${name})
endfunction()

set(PROGRAMS 
	test_constants
	test_math_core
	test_math_calculus
    test_matrix0
    test_matrix1
    test_matrix2
    test_matrix3
    test_matrix4
)

foreach(program ${PROGRAMS})
    add_numlib_test(${program})
	if(WIN32)
		set_tests_properties(${program} PROPERTIES ENVIRONMENT "PATH=${BLAS_BIN};$ENV{PATH}")
	endif()
endforeach()
