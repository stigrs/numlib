cmake_minimum_required(VERSION 3.4.3)
project(liblsode Fortran)

set(
    SRC_FILES
    dlsode.f
)
add_library(
    ${LSODE_LIBRARIES}
    ${SRC_FILES}
)
install(
    TARGETS ${LSODE_LIBRARIES}
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
)

