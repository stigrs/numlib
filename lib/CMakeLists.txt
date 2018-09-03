cmake_minimum_required(VERSION 3.4.3)
project(numlib_lib CXX)

include_directories(${numlib_SOURCE_DIR}/include)
include_directories(${BLAS_INCLUDE_DIRS})
include_directories(${LAPACKE_INCLUDE_DIRS})
link_directories(${numlib_BINARY_DIR})

set(
    SRC_FILES
	calculus.cpp
	linalg.cpp
)
add_library(
    numlib 
    ${SRC_FILES}
)
target_link_libraries(
    numlib 
	${BLAS_LIBRARIES}
	${LAPACK_LIBRARIES}
) 

install(
    TARGETS numlib
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
)
