cmake_minimum_required(VERSION 3.4.3)
project(libquadpack Fortran)

set(
    SRC_FILES
    quadpack.f
)
add_library(
    ${QUADPACK_LIBRARIES}
    ${SRC_FILES}
)

install(
    TARGETS ${QUADPACK_LIBRARIES}
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
)

