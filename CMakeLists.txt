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
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_RELEASE} -std=c++14 -fno-strict-aliasing -Wall -Wcast-align -Wctor-dtor-privacy -Werror -Wextra -Wno-missing-braces -Wno-unknown-attributes -Wnon-virtual-dtor -Wold-style-cast -Woverloaded-virtual -Wpedantic -Wshadow")
endif()
if(MSVC)
    set(CMAKE_CXX_FLAGS_RELEASE "/std:c++14 /W4 /WX /wd4100 /GR /EHsc /MT /O2")
endif()

# Set debugging compiler options.
if (CMAKE_CXX_COMPILER_ID MATCHES "Clang" OR CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    set(CMAKE_CXX_FLAGS_DEBUG "-std=c++14 -g -fno-strict-aliasing -Wall -Wcast-align -Wctor-dtor-privacy -Werror -Wextra -Wno-missing-braces -Wno-unknown-attributes -Wnon-virtual-dtor -Wold-style-cast -Woverloaded-virtual -Wpedantic -Wshadow")
endif()
if(MSVC)
    set(CMAKE_CXX_FLAGS_DEBUG "/std:c++14 /EHsc /W4 /WX /wd4100 /MT /GR /Zi /Od")
endif()

# BLAS and LAPACK is required.
if(DEFINED ENV{MKLROOT})
	list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
	find_package(MKL)
	if(MKL_FOUND)
		if(NOT WIN32)
			if(CMAKE_SIZE_OF_VOID_P EQUAL 8)
				if (CMAKE_CXX_COMPILER_ID MATCHES "Clang" OR CMAKE_CXX_COMPILER_ID MATCHES "GNU")
					set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -m64")
				endif()
			endif()
		endif()
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

list(APPEND GSL_CMAKE_ARGS
    "-DCMAKE_INSTALL_PREFIX=${CMAKE_BINARY_DIR}/external"
    "-DGSL_TEST=OFF"
)

if(GIT_FOUND)
	# Add GSL
	ExternalProject_Add(
		GSL 
		PREFIX ${CMAKE_BINARY_DIR}/GSL 
		GIT_REPOSITORY https://github.com/Microsoft/GSL.git 
		GIT_TAG v2.0.0 
		CMAKE_ARGS ${GSL_CMAKE_ARGS}
		LOG_DOWNLOAD 1
		UPDATE_DISCONNECTED 1
	)
else()
	# Assume GSL is installed in a system directory
	add_custom_target(GSL)
endif()

# Prepare GSL for interface with Numlib.
add_library(gsl_config INTERFACE)
target_include_directories(
	gsl_config INTERFACE
	${CMAKE_BINARY_DIR}/external/include 
)

# Create Numlib interface library.
add_library(numlib_config INTERFACE)
add_dependencies(numlib_config GSL)

# Add definitions.
target_compile_definitions(numlib_config INTERFACE
	GSL_THROW_ON_CONTRACT_VIOLATION
)

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
