// Copyright (c) 2018 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#ifndef NUMLIB_MATH_SPECIAL_FUNCTIONS_H
#define NUMLIB_MATH_SPECIAL_FUNCTIONS_H

#include <numlib/traits.h> 
#include <array>
#include <cmath>
#include <cassert>

namespace Numlib {

// Compute complete elliptic integral of first kind.
template <typename T>
inline Enable_if<Real_type<T>(), T> comp_ellint_1(const T& m)
{
    // Algorithm: Chapter 17.3.34 in Abramovitch and Stegun

    constexpr std::array<T, 5> a = {1.38629436112, 0.09666344259, 0.03590092383,
                                    0.03742563713, 0.01451196212};
    constexpr std::array<T, 5> b = {0.5, 0.12498593597, 0.06880248576,
                                    0.03328355346, 0.00441787012};

    constexpr T zero = T{0};
    constexpr T one = T{1};
    const T m1 = one - m;

    assert(m >= zero && m < one);

    T res = zero;
    for (int i = 0; i < 5; ++i) {
        res += b[i] * std::pow(m1, i);
    }
    res *= std::log(one / m1);
    for (int i = 0; i < 5; ++i) {
        res += a[i] * std::pow(m1, i);
    }
    return res;
}

// Compute complete elliptic integral of second kind.
template <typename T>
inline Enable_if<Real_type<T>(), T> comp_ellint_2(const T& m)
{
    // Algorithm: Chapter 17.3.36 in Abramovitch and Stegun

    constexpr std::array<T, 5> a = {1.0, 0.44325141463, 0.06260601220,
                                    0.04757383546, 0.01736506451};
    constexpr std::array<T, 5> b = {0.0, 0.24998368310, 0.09200180037,
                                    0.04069697526, 0.00526449639};

    constexpr T zero = T{0};
    constexpr T one = T{1};
    const T m1 = one - m;

    assert(m >= zero && m < one);

    T res = zero;
    for (int i = 0; i < 5; ++i) {
        res += b[i] * std::pow(m1, i);
    }
    res *= std::log(one / m1);
    for (int i = 0; i < 5; ++i) {
        res += a[i] * std::pow(m1, i);
    }
    return res;
}

} // namespace Numlib

#endif // NUMLIB_MATH_SPECIAL_FUNCTIONS_H

