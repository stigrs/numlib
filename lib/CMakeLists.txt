cmake_minimum_required(VERSION 3.4.3)
project(libnumlib CXX)

include_directories(${BLAS_INCLUDE_DIRS})
include_directories(${LAPACKE_INCLUDE_DIRS})
link_directories(${numlib_BINARY_DIR})
if(MKL_FOUND)
    link_directories(${MKL_LIBRARY_DIRS})
endif()

set(
    SRC_FILES
	calculus.cpp
    geometry.cpp
	linalg.cpp
)
add_library(
    numlib 
    ${SRC_FILES}
)
add_dependencies(numlib numlib_config)
target_link_libraries(
    numlib 
	${BLAS_LIBRARIES}
    ${LAPACK_LIBRARIES}
    numlib_config
) 

install(
    TARGETS numlib
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
)
