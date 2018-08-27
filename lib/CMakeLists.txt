cmake_minimum_required(VERSION 3.4.3)
project(numlib_lib CXX)

include_directories(${numlib_SOURCE_DIR}/include)
if(NOT APPLE)
	include_directories(${BLAS_INCLUDE_DIR})
	include_directories(${LAPACKE_INCLUDE_DIR})
endif()

link_directories(${numlib_BINARY_DIR})

set(
    SRC_FILES
	calculus.cpp
)
add_library(
    numlib 
    ${SRC_FILES}
)
target_link_libraries(
    numlib 
) 

install(
    TARGETS numlib
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
)
