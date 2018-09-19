// Copyright (c) 2018 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#ifndef NUMLIB_MATH_CORE_H
#define NUMLIB_MATH_CORE_H

#include <numlib/traits.h>
#include <cstdlib>
#include <cmath>

namespace Numlib {

//------------------------------------------------------------------------------
//
// Provides core mathematical functions:

// Check if integer type is even.
template <typename T>
inline Enable_if<Integer_type<T>(), bool> even(const T& n)
{
    return n % 2 ? false : true;
}

// Check if integer type is odd.
template <typename T>
inline Enable_if<Integer_type<T>(), bool> odd(const T& n)
{
    return n % 2 ? true : false;
}

// Compute Kronecker delta.
template <typename T>
inline Enable_if<Integer_type<T>(), T> krond(const T& i, const T& j)
{
    return i == j;
}

// Round double to nearest integer type.
template <typename T>
inline Enable_if<Integer_type<T>(), T> round(double x)
{
    return static_cast<T>((x > 0) ? x + 0.5 : x - 0.5);
}

// Sign transfer function.
template <typename T>
inline T sign(const T& x, const T& y)
{
    return (y >= 0) ? std::abs(x) : -std::abs(x);
}

// Raise argument to the power of 2.
template <typename T>
inline T sqr(const T& x)
{
    return x * x;
}

} // namespace Numlib

#endif // NUMLIB_MATH_CORE_H
