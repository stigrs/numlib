cmake_minimum_required(VERSION 3.4.3)
project(numlib CXX)

include(ExternalProject)
find_package(Git)

# By default, numlib is install in the HOME directory.
if(WIN32)
    set(CMAKE_INSTALL_PREFIX "C:/Users/$ENV{USERNAME}")
else()
    set(CMAKE_INSTALL_PREFIX $ENV{HOME})
endif()

# Enforce C++14 standard.
set(CMAKE_CXX_STANDARD 14)

# Export compile commands.
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(WIN32)
    option(BUILD_SHARED_LIBS "Compile shared libraries." OFF)
else()
    option(BUILD_SHARED_LIBS "Compile shared libraries." ON)
endif()

option(BUILD_TESTS "Build tests." ON)
option(BUILD_BENCH "Build benchmarks." OFF)
if(WIN32)
	set(BUILD_BENCH "Build benchmarks." ON)
endif()

# Set default MSVC compiler options to avoid D9025 error.
if(MSVC)
    set(CMAKE_CXX_FLAGS "/DWIN32 /D_WINDOWS")
endif()

# Set release compiler options.
if (CMAKE_CXX_COMPILER_ID MATCHES "Clang" OR CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_RELEASE} -Wpedantic -Wall -Wshadow -Wextra -Wno-missing-braces -Werror -Wno-error=strict-overflow")
endif()
if(MSVC)
    set(CMAKE_CXX_FLAGS_RELEASE "/DNDEBUG /std:c++14 /W4 /WX /wd4100 /GR /EHsc /MT /O2")
endif()

# Set debugging compiler options.
if (CMAKE_CXX_COMPILER_ID MATCHES "Clang" OR CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wpedantic -Wall -Wshadow -Wextra -Wno-missing-braces -Werror -Wno-error=strict-overflow")
endif()
if(MSVC)
    set(CMAKE_CXX_FLAGS_DEBUG "/std:c++14 /EHsc /W4 /WX /wd4100 /MT /GR /Zi /Od")
endif()

# BLAS and LAPACK is required.
if(DEFINED ENV{MKLROOT})
	list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
	find_package(MKL)
	if(MKL_FOUND)
        #		if(NOT WIN32)
        #	if(CMAKE_SIZEOF_VOID_P EQUAL 8)
        #		if (CMAKE_CXX_COMPILER_ID MATCHES "Clang" OR CMAKE_CXX_COMPILER_ID MATCHES "GNU")
        #			set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -m64")
        #		endif()
        #	endif()
        #endif()
		add_definitions(${MKL_DEFINITIONS})
		set(BLAS_INCLUDE_DIRS ${MKL_INCLUDE_DIRS})
		set(BLAS_LIBRARIES ${MKL_LIBRARIES})
		set(LAPACK_LIBRARIES "")
	endif()
else()
	set(BLA_VENDOR "OpenBLAS")
	if(WIN32)
		set(BLAS_DIR C:/local/OpenBLAS.0.2.14.1/lib/native)
		if(${CMAKE_CL_64} EQUAL 1)
			find_file(BLAS_LIBRARIES libopenblas.dll.a HINTS ${BLAS_DIR}/lib/x64)
			set(BLAS_BIN "${BLAS_DIR}/bin/x64")
		else()
			find_file(BLAS_LIBRARIES libopenblas.dll.a HINTS ${BLAS_DIR}/lib/win32)
			set(BLAS_BIN "${BLAS_DIR}/bin/win32")
		endif()
	elseif(APPLE)
		set(BLAS_DIR /usr/local/opt/openblas)
		find_file(BLAS_LIBRARIES libopenblas.dylib HINTS ${BLAS_DIR}/lib)
	else()
		# Workaround: LAPACKE is not a part of OpenBLAS on Ubuntu
		set(BLAS_DIR /usr)
		find_package(BLAS)
		find_file(LAPACK_LIBRARIES liblapacke.so HINTS ${BLAS_DIR}/lib ${BLAS_DIR}/lib/x86_64-linux-gnu) 
		find_path(BLAS_INCLUDE_DIRS cblas.h HINTS ${BLAS_DIR}/include /usr/include /usr/local/include)
		find_path(LAPACKE_INCLUDE_DIRS lapacke.h HINTS ${BLAS_DIR}/include /usr/include /usr/local/include)
	endif()
	find_path(BLAS_INCLUDE_DIRS cblas.h HINTS ${BLAS_DIR}/include)
	find_path(LAPACKE_INCLUDE_DIRS lapacke.h HINTS ${BLAS_DIR}/include)
endif()

# Create Numlib interface library.
add_library(numlib_config INTERFACE)
target_include_directories(numlib_config INTERFACE
    $<BUILD_INTERFACE:
        ${CMAKE_CURRENT_SOURCE_DIR}/include
    >
)

install(
    DIRECTORY include/numlib
    DESTINATION include
    FILES_MATCHING PATTERN "*.h"
)

add_subdirectory(lib)

if(BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()

if(BUILD_BENCH)
    # Armadillo library is required.
    if(WIN32)
		set(ARMA_DIR C:/local/armadillo-code.7.900.1/build/native)
        find_file(ARMADILLO_LIBRARIES armadillo.lib HINTS ${ARMA_DIR}/lib)
        find_path(ARMADILLO_INCLUDE_DIRS armadillo HINTS ${ARMA_DIR}/include)
        add_definitions(-DARMA_NO_DEBUG)
    else()
        find_package(Armadillo 7.900.1 REQUIRED)
        add_definitions(-DARMA_NO_DEBUG)
    endif()
	add_subdirectory(bench)
endif()
