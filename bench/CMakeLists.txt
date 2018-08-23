cmake_minimum_required(VERSION 3.4.3)
project(numlib_bench CXX)

include_directories(${numlib_SOURCE_DIR}/include)
include_directories(${ARMADILLO_INCLUDE_DIRS})
if(NOT APPLE)
	include_directories(${BLAS_INCLUDE_DIR})
	include_directories(${LAPACKE_INCLUDE_DIR})
endif()

set(
    PROGRAMS 
    bench_axpy 
	bench_mm_mul
)

foreach(program ${PROGRAMS})
    add_executable(${program} ${program}.cpp)
	if(WIN32)
		target_link_libraries (
			${program} 
			${BLAS_LIBRARIES}
			${ARMADILLO_LIBRARIES} 
		) 
	else()
		target_link_libraries (
			${program} 
			${BLAS_LIBRARIES}
		) 
	endif()
endforeach()
