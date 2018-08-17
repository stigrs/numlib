cmake_minimum_required(VERSION 3.4.3)
project(numlib_tests CXX)

include_directories(${numlib_SOURCE_DIR}/include)
include_directories(${GSL_INCLUDE_DIR})
include_directories(${ARMADILLO_INCLUDE_DIRS})
include_directories(${MKL_INCLUDE_DIR})
link_directories(${MKL_LIB_DIR})
if(APPLE)
    link_directories(${TBB_LIB_DIR})
endif()

# Copy *.inp files to build directory.
file(GLOB INP_FILES RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}/*.inp")

foreach(inp_file_name ${INP_FILES})
    file(COPY "${inp_file_name}" DESTINATION "${CMAKE_BINARY_DIR}/tests")
endforeach()

# Prepare "Catch" library for other executables.
set(CATCH_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR})
add_library(Catch INTERFACE)
target_include_directories(Catch INTERFACE ${CATCH_INCLUDE_DIR})

# Create the main executable for each test.
add_library(test_catch STATIC test.cpp)
target_link_libraries(
    test_catch
    Catch
)

function(add_numlib_test name)
    add_executable(${name} ${name}.cpp)
    if(WIN32)
        target_link_libraries (
            ${program} 
            mkl_intel_c.lib 
            mkl_sequential.lib 
            mkl_core.lib 
            test_catch
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
            test_catch
        )
    endif()
    add_dependencies(${name} Catch)
    add_test(NAME ${name} COMMAND ${name})
endfunction()

set(PROGRAMS 
    test_matrix1
    test_matrix2
    test_matrix3
)

foreach(program ${PROGRAMS})
    add_numlib_test(${program})
endforeach()
