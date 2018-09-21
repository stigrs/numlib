cmake_minimum_required(VERSION 3.4.3)
project(numlib_examples CXX)

include_directories(${BLAS_INCLUDE_DIRS})
include_directories(${LAPACKE_INCLUDE_DIRS})
if(MKL_FOUND)
	link_directories(${MKL_LIBRARY_DIRS})
endif()

set(
    PROGRAMS 
    arithmetic
    band_matrix
    constants
    convolution
    derivation
    eigensolvers
    integration
    linear_algebra
    matrix_decomposition 
    matrix_inversion 
    matrix_properties 
    packed_matrix
    solve_linear_equations
    sparse_eigensolver
    sparse_linsolve
    statistics
	sub_matrix_views
	subscripting_and_slicing
)

foreach(program ${PROGRAMS})
    add_executable(${program} ${program}.cpp)
	target_link_libraries (
		${program} 
		numlib
		${BLAS_LIBRARIES}
		${LAPACK_LIBRARIES}
		numlib_config
	) 
endforeach()
