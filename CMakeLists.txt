cmake_minimum_required(VERSION 3.4.3)
project(numlib CXX)

# By default, numlib is install in the HOME directory.
if(WIN32)
    set(CMAKE_INSTALL_PREFIX "C:/Users/$ENV{USERNAME}")
else()
    set(CMAKE_INSTALL_PREFIX $ENV{HOME})
endif()

# Export compile commands.
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Enforce C++14 standard (needed for GSL).
set(CMAKE_CXX_STANDARD 14)

if(WIN32)
    option(BUILD_SHARED_LIBS "Compile shared libraries." OFF)
else()
    option(BUILD_SHARED_LIBS "Compile shared libraries." ON)
endif()

option(BUILD_TESTS "Build tests." ON)

# Set default MSVC compiler options to avoid D9025 error.
if(MSVC)
    set(CMAKE_CXX_FLAGS "/DWIN32 /D_WINDOWS")
endif()

# Set release compiler options.
if (CMAKE_CXX_COMPILER_ID MATCHES "Clang" OR CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -Wpedantic -Wall -Wshadow -Wextra -Wno-missing-braces -Werror -Wno-error=strict-overflow")
endif()
if (CMAKE_CXX_COMPILER_ID MATCHES "Intel")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -std=c++14 -Wall -Wshadow -Wextra -Wno-missing-braces -Werror -Wno-error=strict-overflow")
endif()
if(APPLE)
    set(CMAKE_EXE_LINKER_FLAGS_RELEASE "${CMAKE_EXE_LINKER_FLAGS_RELEASE} -framework Accelerate")
endif()
if(MSVC)
    set(CMAKE_CXX_FLAGS_RELEASE "/std:c++14 /W4 /WX /wd4100 /GR /EHsc /MT /O2")
endif()

# Set debugging compiler options.
if (CMAKE_CXX_COMPILER_ID MATCHES "Clang" OR CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -pg -Wpedantic -Wall -Wshadow -Wextra -Wno-missing-braces -Werror -Wno-error=strict-overflow")
    set(CMAKE_EXE_LINKER_FLAGS_DEBUG "${CMAKE_EXE_LINKER_FLAGS_DEBUG} -pg")
endif()
if (CMAKE_CXX_COMPILER_ID MATCHES "Intel") 
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -std=c++14 -pg -Wall -Wshadow -Wextra -Wno-missing-braces -Werror -Wno-error=strict-overflow")
    set(CMAKE_EXE_LINKER_FLAGS_DEBUG "${CMAKE_EXE_LINKER_FLAGS_DEBUG} -pg")
endif()
if(MSVC)
    set(CMAKE_CXX_FLAGS_DEBUG "/std:c++14 /EHsc /W4 /WX /wd4100 /MT /GR /Zi /Od")
endif()

# Set definitions.
set(numlib_DEFINITIONS_DEBUG -DGSL_THROW_ON_CONTRACT_VIOLATION)
set(numlib_DEFINITIONS_RELEASE -DNDEBUG -DGSL_THROW_ON_CONTRACT_VIOLATION)
add_compile_options(
    "$<$<CONFIG:Debug>:${numlib_DEFINITIONS_DEBUG}>"
    "$<$<CONFIG:Release>:${numlib_DEFINITIONS_RELEASE}>" 
)

# Boost library is required.
#if(WIN32)
#    set(Boost_USE_STATIC_LIBS ON)
#endif()
#find_package(Boost 1.65.0 REQUIRED)

# GSL library is required.
#find_path(GSL_INCLUDE_DIR gsl HINTS $ENV{HOME}/include /usr/include /usr/local/include)

# Intel MKL is required.
#find_path(MKL_INCLUDE_DIR mkl.h HINTS $ENV{MKL_INCLUDE_DIR} $ENV{MKLROOT}/include /opt/intel/mkl/include)
#if (WIN32)
#    find_path(MKL_LIB_DIR mkl_core.lib HINTS $ENV{MKL_LIB_DIR})
#else()
#    find_path(MKL_LIB_DIR libmkl_core.a HINTS $ENV{MKL_LIB_DIR} $ENV{MKLROOT}/lib/intel64 /opt/intel/mkl/lib /opt/intel/mkl/lib/intel64)
#    if(APPLE)
#        find_path(TBB_LIB_DIR libtbb.dylib HINTS /opt/intel/tbb/lib)
#    endif()
#endif()

add_library(numlib_h INTERFACE)

target_include_directories(numlib_h INTERFACE
    $<BUILD_INTERFACE:
        ${CMAKE_CURRENT_SOURCE_DIR}/include
        ${GSL_INCLUDE_DIR}
    >
)

install(
    DIRECTORY include/numlib
    DESTINATION include
    FILES_MATCHING PATTERN "*.h"
)

#add_subdirectory(lib)

if(BUILD_TESTS)
    # Armadillo library is required.
#    if(WIN32)
#        find_path(ARMADILLO_INCLUDE_DIRS armadillo HINTS $ENV{HOME}/include)
#        add_definitions(-DARMA_USE_BLAS)
#        add_definitions(-DARMA_USE_LAPACK)
#        add_definitions(-DARMA_NO_DEBUG)
#    else()
#        find_package(Armadillo 7.950.1 REQUIRED)
#        add_definitions(-DARMA_NO_DEBUG)
#    endif()
    enable_testing()
    add_subdirectory(tests)
#   add_subdirectory(bench)
endif()
