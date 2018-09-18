cmake_minimum_required(VERSION 3.4.3)
project(numlib_examples CXX)

include_directories(${BLAS_INCLUDE_DIRS})
include_directories(${LAPACKE_INCLUDE_DIRS})
if(MKL_FOUND)
	link_directories(${MKL_LIBRARY_DIRS})
endif()

set(
    PROGRAMS 
    basic_matrix_uses 
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
