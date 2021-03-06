# numlib [![Build Status](https://travis-ci.org/stigrs/numlib.svg?branch=master)](https://travis-ci.org/stigrs/numlib)[![Build status](https://ci.appveyor.com/api/projects/status/github/stigrs/numlib?svg=true)](https://ci.appveyor.com/project/stigrs/numlib)

Numlib provides a C++ library for linear algebra and scientific computing.
BLAS and LAPACK are used for fast numerical performance. Currently, OpenBLAS
and Intel MKL are supported.

## Features

* N-dimensional dense matrices using Stroustrup's matrix design
  (row-major storage order)
* Band matrices (column-major storage order)
* Packed matrices (row-major storage order)
* Sparse vectors and matrices (CSR3 storage format)
* Selected mathematical functions not provided by the STL
* Faddeeva package (w, erf, erfc, erfcx, erfi, Dawson)
* Numerical derivation and integration methods (including QAGS and QAGI from 
  CQUADPACK)
* Solvers for ordinary differential equations (DOPRI5 and LSODA)
* Linear algebra methods
* Transformations between rotation matrix, Euler angles and quaternions
* Vector convolution
* Forward and inverse discrete Fourier transform (FFT and IFFT)
* Generation of mesh grids
* Statistical methods
* Mathematical constants, metric prefixes, physical constants, and
  conversion factors

_Note: Some features are only available if Intel MKL is used._

## Code of Conduct

This project has adopted the [Covenant Code of Conduct](CODE_OF_CONDUCT.md).

## Licensing

Numlib is released under the [MIT](LICENSE) license.

## Usage of Third Party Libraries

This project makes use of code from the following third-party libraries:

* [Catch2](https://github.com/catchorg/Catch2)
* [origin](http://code.google.com/p/origin)
* [Faddeeva package](http://ab-initio.mit.edu/wiki/index.php/Faddeeva_Package) 
* [CQUADPACK](https://github.com/ESSS/cquadpack.git) 
* [LSODA](https://github.com/lh3/misc.git)

Please see the [ThirdPartyNotices.txt](ThirdPartyNotices.txt) file for details
regarding the licensing of these libraries.

## Quick Start

### Requirements

* [CMake](https://cmake.org) 3.4.3
* [Boost](https://www.boost.org) 1.54.0
* [OpenBLAS](https://www.openblas.net/) 0.3.3 (Intel MKL is recommended)
* [Armadillo](http://arma.sourceforge.net) 7.900.1 (for benchmarking)

### Supported Platforms

The test suite that exercises Numlib has been built and passes successfully
on the following platforms:
* GNU/Linux using GCC 6.4.0, 7.3.0
* GNU/Linux using Clang 6.0
* OS X High Sierra (10.13) using Apple Xcode 9.4, 10.0 
* macOS Mojave (10.14) using Apple Xcode 11.3 
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
