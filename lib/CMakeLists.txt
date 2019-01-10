cmake_minimum_required(VERSION 3.4.3)
project(libnum CXX)

include_directories(${numlib_SOURCE_DIR}/include)
include_directories(${Boost_INCLUDE_DIRS})
include_directories(${BLAS_INCLUDE_DIRS})
include_directories(${LAPACKE_INCLUDE_DIRS})
link_directories(${numlib_BINARY_DIR})
if(MKL_FOUND)
    link_directories(${MKL_LIBRARY_DIRS})
endif()
if(ENABLE_QUADPACK)
    link_directories(${CMAKE_BINARY_DIR}/quadpack)
endif()

set(
    SRC_FILES
    calculus.cpp
    euler.cpp
    faddeeva.cpp
    geometry.cpp
    grid.cpp
    linalg.cpp
    quaternion.cpp
    signal.cpp
    statistics.cpp
)
add_library(
    num 
    ${SRC_FILES}
)
target_link_libraries(
    num 
    ${BLAS_LIBRARIES}
    ${QUADPACK_LIBRARIES}
) 

install(
    TARGETS num
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
)

