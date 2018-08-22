cmake_minimum_required(VERSION 3.4.3)
project(numlib_bench CXX)

include_directories(${numlib_SOURCE_DIR}/include)
include_directories(${ARMADILLO_INCLUDE_DIRS})
include_directories(${BLAS_INCLUDE_DIR})
include_directories(${LAPACKE_INCLUDE_DIR})

set(
    PROGRAMS 
    bench_axpy 
	bench_mm_mul
)

foreach(program ${PROGRAMS})
    add_executable(${program} ${program}.cpp)
    target_link_libraries (
        ${program} 
		${BLAS_LIBRARIES}
        ${ARMADILLO_LIBRARIES} 
    ) 
endforeach()
