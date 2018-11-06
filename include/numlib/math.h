// Copyright (c) 2018 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#ifndef NUMLIB_MATH_H
#define NUMLIB_MATH_H

#include <stdexcept>
#include <string>

#ifdef _MSC_VER // ugly hack since Visual Studio does not support C99 _Complex
#ifndef USE_MKL
#include <complex>
#define lapack_complex_float std::complex<float>
#define lapack_complex_double std::complex<double>
#endif
#endif

// Math error.
struct Math_error : std::runtime_error {
    Math_error(const std::string& s) : std::runtime_error(s) {}
};

// Note: Do not change this order.
#include <numlib/math_impl/core.h>
#include <numlib/math_impl/calculus.h>
#include <numlib/math_impl/linalg.h>
#include <numlib/math_impl/geometry.h>
#include <numlib/math_impl/euler.h>
#include <numlib/math_impl/quaternion.h>
#include <numlib/math_impl/signal.h>
#include <numlib/math_impl/special_functions.h>
#include <numlib/math_impl/statistics.h>

#endif // NUMLIB_MATH_H
