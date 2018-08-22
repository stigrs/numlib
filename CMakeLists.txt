cmake_minimum_required(VERSION 3.7.2)
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
set(numlib_DEFINITIONS_DEBUG)
set(numlib_DEFINITIONS_RELEASE -DNDEBUG)
add_compile_options(
    "$<$<CONFIG:Debug>:${numlib_DEFINITIONS_DEBUG}>"
    "$<$<CONFIG:Release>:${numlib_DEFINITIONS_RELEASE}>" 
)

# BLAS and LAPACK is required.
set(BLA_VENDOR "OpenBLAS")
if(APPLE)
    set(BLA_VENDOR "Apple")
endif()
if(WIN32)
    set(BLAS_DIR C:/local/OpenBLAS.0.2.14.1/lib/native)
	if(${CMAKE_CL_64} EQUAL 1)
		find_file(BLAS_LIBRARIES libopenblas.dll.a HINTS ${BLAS_DIR}/lib/x64)
		set(BLAS_BIN "${BLAS_DIR}/bin/x64")
	else()
		find_file(BLAS_LIBRARIES libopenblas.dll.a HINTS ${BLAS_DIR}/lib/win32)
		set(BLAS_BIN "${BLAS_DIR}/bin/win32")
	endif()
    find_path(BLAS_INCLUDE_DIR cblas.h HINTS ${BLAS_DIR}/include)
    find_path(LAPACKE_INCLUDE_DIR lapacke.h HINTS ${BLAS_DIR}/include)
else()
	find_package(BLAS REQUIRED HINTS})
    find_path(BLAS_INCLUDE_DIR cblas.h HINTS ${BLAS_DIR}/include /usr/include /usr/local/include)
    find_path(LAPACKE_INCLUDE_DIR lapacke.h HINTS ${BLAS_DIR}/include /usr/include /usr/local/include)
endif()

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
    if(WIN32)
		set(ARMA_DIR C:/local/armadillo-code.7.900.1/build/native)
        find_file(ARMADILLO_LIBRARIES armadillo.lib HINTS ${ARMA_DIR}/lib)
        find_path(ARMADILLO_INCLUDE_DIRS armadillo HINTS ${ARMA_DIR}/include)
        add_definitions(-DARMA_NO_DEBUG)
    endif()
    enable_testing()
    add_subdirectory(tests)
	if(WIN32)
		add_subdirectory(bench)
	endif()
endif()
