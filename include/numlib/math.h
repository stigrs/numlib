// Copyright (c) 2018 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#ifndef NUMLIB_MATH_H
#define NUMLIB_MATH_H

#include <stdexcept>
#include <string>

// Math error.
struct Math_error : std::runtime_error {
    Math_error(const std::string& s) : std::runtime_error(s) {}
};

#include <numlib/math_impl/core.h>
#include <numlib/math_impl/calculus.h>
#include <numlib/math_impl/linalg.h>

#endif // NUMLIB_MATH_H
