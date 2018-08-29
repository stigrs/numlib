cmake_minimum_required(VERSION 3.4.3)
project(numlib_bench CXX)

include_directories(${numlib_SOURCE_DIR}/include)
if(WIN32)
    include_directories(${ARMADILLO_INCLUDE_DIRS})
endif()
if(NOT APPLE)
	include_directories(${BLAS_INCLUDE_DIR})
	include_directories(${LAPACKE_INCLUDE_DIR})
endif()

set(
    PROGRAMS 
    bench_axpy 
    bench_dot 
	bench_mm_mul
	bench_mv_mul
	bench_norm
)

foreach(program ${PROGRAMS})
    add_executable(${program} ${program}.cpp)
	target_link_libraries (
		${program} 
		${BLAS_LIBRARIES}
		${LAPACK_LIBRARIES}
		${ARMADILLO_LIBRARIES} 
	) 
endforeach()
