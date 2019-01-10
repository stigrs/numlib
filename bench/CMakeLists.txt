cmake_minimum_required(VERSION 3.4.3)
project(numlib_bench CXX)

if(WIN32)
    include_directories(${ARMADILLO_INCLUDE_DIRS})
endif()
include_directories(${numlib_SOURCE_DIR}/include)
include_directories(${Boost_INCLUDE_DIRS})
include_directories(${BLAS_INCLUDE_DIRS})
include_directories(${LAPACKE_INCLUDE_DIRS})
if(MKL_FOUND)
    link_directories(${MKL_LIBRARY_DIRS})
endif()

set(
    PROGRAMS 
    bench_axpy 
    bench_dot 
    bench_eigs
    bench_mm_mul
    bench_mv_mul
    bench_norm
    bench_transpose
)

foreach(program ${PROGRAMS})
    add_executable(${program} ${program}.cpp)
	target_link_libraries (
	    ${program} 
	    num
	    ${BLAS_LIBRARIES}
	    ${LAPACK_LIBRARIES}
	    ${ARMADILLO_LIBRARIES} 
	) 
endforeach()
