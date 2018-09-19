cmake_minimum_required(VERSION 3.4.3)
project(numlib_tests CXX)

include_directories(${BLAS_INCLUDE_DIRS})
include_directories(${LAPACKE_INCLUDE_DIRS})
if(MKL_FOUND)
	link_directories(${MKL_LIBRARY_DIRS})
endif()

list(APPEND CATCH_CMAKE_ARGS
    "-DCMAKE_INSTALL_PREFIX=${CMAKE_BINARY_DIR}/external"
    "-DCATCH_BUILD_TESTING=OFF"
)

if(GIT_FOUND)
    # Add catch
    ExternalProject_Add(
        catch2
        PREFIX ${CMAKE_BINARY_DIR}/catch
        GIT_REPOSITORY https://github.com/catchorg/Catch2.git
        GIT_TAG v2.3.0
        CMAKE_ARGS ${CATCH_CMAKE_ARGS}
        LOG_DOWNLOAD 1
        UPDATE_DISCONNECTED 1
    )
else()
    # Assume catch is installed in a system directory
    add_custom_target(catch2)
endif()

# Prepare "Catch" library for other executables.
add_library(Catch2 INTERFACE)
target_include_directories(
    Catch2 INTERFACE 
    ${CMAKE_BINARY_DIR}/external/include
)

# Create the main executable for each test.
add_library(test_catch2 STATIC test.cpp)
target_link_libraries(
    test_catch2
    Catch2
)

add_dependencies(test_catch2 catch2)

function(add_numlib_test name)
    add_executable(${name} ${name}.cpp)
    target_link_libraries (
        ${name} 
		numlib
		${BLAS_LIBRARIES}
		${LAPACK_LIBRARIES}
        test_catch2
        numlib_config
    ) 
    add_dependencies(${name} Catch2)
    add_test(NAME ${name} COMMAND ${name})
endfunction()

set(PROGRAMS 
	test_band_matrix
	test_constants
	test_math_core
	test_math_calculus
    test_math_euler
	test_math_geometry
	test_math_linalg
	test_math_quaternion
	test_math_signal
	test_math_statistics
    test_matrix0
    test_matrix1
    test_matrix2
    test_matrix3
	test_matrix4
	test_packed_matrix
	test_sparse_matrix
	test_sparse_vector
)

foreach(program ${PROGRAMS})
    add_numlib_test(${program})
	if(WIN32)
		if(NOT MKL_FOUND)
			set_tests_properties(${program} PROPERTIES ENVIRONMENT "PATH=${BLAS_BIN};$ENV{PATH}")
		endif()
	endif()
endforeach()
