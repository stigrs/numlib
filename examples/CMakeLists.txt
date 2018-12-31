cmake_minimum_required(VERSION 3.4.3)
project(numlib_examples CXX)

include_directories(${numlib_SOURCE_DIR}/include)
include_directories(${BLAS_INCLUDE_DIRS})
include_directories(${LAPACKE_INCLUDE_DIRS})
if(MKL_FOUND)
    link_directories(${MKL_LIBRARY_DIRS})
endif()
if(BUILD_LSODE)
    link_directories(${CMAKE_BINARY_DIR}/liblsode)
endif()

set(
    PROGRAMS 
    arithmetic
    band_matrix
    constants
    convolution
    derivation
    eigensolvers
    fourier_transform
    integration
    linear_algebra
    linear_least_squares
    lsode_lorenz
    lsode_stiff
    matrix_decomposition 
    matrix_inversion 
    matrix_properties 
    packed_matrix
    solve_linear_equations
    sparse_eigensolver
    sparse_linsolve
    special_functions
    statistics
    sub_matrix_views
    subscripting_and_slicing
)

foreach(program ${PROGRAMS})
    add_executable(${program} ${program}.cpp)
	target_link_libraries (
	    ${program} 
	    num
	    ${BLAS_LIBRARIES}
	    ${LAPACK_LIBRARIES}
        ${LSODE_LIBRARIES}
	) 
endforeach()
