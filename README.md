# numlib [![Build Status](https://travis-ci.org/stigrs/numlib.svg?branch=master)](https://travis-ci.org/stigrs/numlib)[![Build status](https://ci.appveyor.com/api/projects/status/github/stigrs/numlib?svg=true)](https://ci.appveyor.com/project/stigrs/numlib)

Numlib provides a C++ library for linear algebra and scientific computing. 
BLAS and LAPACK are used for fast numerical performance. 

## Features

* N-dimensional dense matrices using Stroustrup's matrix design (see TC++PL4) 
* A few basic mathematical functions not provided by the STL
* Basic numerical derivation and integration methods
* Mathematical constants, metric prefixes, physical constants, and 
  conversion factors

## Code of Conduct

This project has adopted the [Covenant Code of Conduct](CODE_OF_CONDUCT.md).

## Licensing

Numlib is released under the [MIT](LICENSE) license.

## Usage of Third Party Libraries

This project makes use of the [Catch2](https://github.com/catchorg/Catch2) 
testing library and code from [origin](http://code.google.com/p/origin). 
Please see the [ThirdPartyNotices.txt](ThirdPartyNotices.txt) file for details 
regarding the licensing of Catch2 and origin.

## Quick Start 

### Requirements

* [CMake](https://cmake.org) 3.4.3
* [OpenBLAS](https://www.openblas.net/) 0.2.14.1 (Linux and Windows)
* [Armadillo](http://arma.sourceforge.net) 7.900.1 (for benchmarking)

### Supported Platforms

The test suite that exercises Numlib has been built and passes successfully 
on the following platforms:
* GNU/Linux using GCC 5.5.0, 6.4.0, 7.3.0
* GNU/Linux using Clang 3.6, 3.7, 3.8, 3.9, 4.0, 5.0
* OS X El Capitan (10.12) using Apple LLVM 8.3.0
* OS X High Sierra (10.13) using Apple LLVM 9.1, 9.4
* Windows using Visual Studio 2017 (x86 and x64)

### Obtaining the Source Code

The source code can be obtained from

        git clone git@github.com:stigrs/numlib.git

### Building the Software

These steps assumes that the source code of this repository has been cloned
into a directory called `numlib`.

1. Create a directory to contain the build outputs:

        cd numlib
        mkdir build
        cd build

2. Configure CMake to use the compiler of your choice (you can see a list by
   running `cmake --help`):

        cmake -G "Visual Studio 15 2017" ..

3. Build the software (in this case in the Release configuration):

        cmake --build . --config Release

4. Run the test suite:

        ctest -C Release

5. Install the software:

        cmake --build . --config Release --target install

   All tests should pass, indicating that your platform is fully supported. 

6. Benchmarks can be built by setting the option BUILD_BENCH to ON. Please 
   make sure BLAS run on the same number of threads in Armadillo and Numlib 
   before comparing the benchmark results. If OpenBLAS is used, this can be
   controlled by setting the OPENBLAS_NUM_THREADS environmental variable.
