cmake_minimum_required(VERSION 3.4.3)
project(libodepack Fortran)

set(
    SRC_FILES
    dlsode.f
    rkf45.f
)
add_library(
    ${ODEPACK_LIBRARIES}
    ${SRC_FILES}
)
install(
    TARGETS ${ODEPACK_LIBRARIES}
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
)

