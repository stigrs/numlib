cmake_minimum_required(VERSION 3.4.3)
project(numlib_bench CXX)

include_directories(${numlib_SOURCE_DIR}/include)
include_directories(${ARMADILLO_INCLUDE_DIRS})
include_directories(${MKL_INCLUDE_DIR})
link_directories(${MKL_LIB_DIR})
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
            mkl_intel_c.lib 
            mkl_sequential.lib 
            mkl_core.lib
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
            test_catch
        )
    else()
        target_link_libraries (
            ${program} 
            mkl_intel_lp64
            mkl_gnu_thread
            mkl_core
            gomp
            pthread
            m 
            dl 
            ${ARMADILLO_LIBRARIES} 
        )
    endif()
endforeach()
