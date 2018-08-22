cmake_minimum_required(VERSION 3.4.3)
project(numlib_bench CXX)

include_directories(${numlib_SOURCE_DIR}/include)
include_directories(${ARMADILLO_INCLUDE_DIRS})
include_directories(${BLAS_INCLUDE_DIR})
include_directories(${LAPACKE_INCLUDE_DIR})
if(APPLE)
    link_directories(${TBB_LIB_DIR})
endif()

set(
    PROGRAMS 
    bench_axpy 
)

foreach(program ${PROGRAMS})
    add_executable(${program} ${program}.cpp)
    if(WIN32)
        target_link_libraries (
            ${program} 
			${BLAS_LIBRARIES}
        ) 
    elseif(APPLE)
        target_link_libraries (
            ${program} 
            mkl_intel_lp64
            mkl_tbb_thread
            mkl_core
            tbb
            stdc++
            pthread
            m 
            ldl 
            ${ARMADILLO_LIBRARIES}
        )
    else()
        target_link_libraries (
            ${program} 
			${BLAS_LIBRARIES}
            ${ARMADILLO_LIBRARIES} 
        )
    endif()
endforeach()
