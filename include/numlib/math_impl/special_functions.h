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
inline Enable_if<Real_type<T>(), T> comp_ellint_1(const T& k)
{
    // Algorithm:
    // ----------
    // Cody, W. J. Chebyshev approximations for the complete elliptic
    // integrals K and E. Mathematics of Computation, 1965, vol. 19,
    // pp. 105-112.
    //
    // Coefficients are taken from Table II.

    assert(k >= 0.0 && k < 1.0);
    const T eta = 1.0 - k * k;

    constexpr int n = 8;
    constexpr std::array<T, n> a = {
        9.6573590797589018e-2, 3.0885573486752694e-2, 1.4978988178704629e-2,
        9.6587579861753112e-3, 1.1208918554644092e-2, 1.3855601247156560e-2,
        6.6905509906897936e-3, 6.4998443329390180e-4};
    constexpr std::array<T, n> b = {
        1.2499999994117923e-1, 7.0312426464627361e-2, 4.8818058565403952e-2,
        3.7068398934155422e-2, 2.7189861116788250e-2, 1.4105380776158048e-2,
        3.1831309927862886e-3, 1.5049181783601883e-4};

    T res = 0.5;
    for (int i = 0; i < n; ++i) {
        res += b[i] * std::pow(eta, i + 1);
    }
    res *= std::log(1.0 / eta);
    for (int i = 0; i < n; ++i) {
        res += a[i] * std::pow(eta, i + 1);
    }
    return res + std::log(4.0);
}

// Compute complete elliptic integral of second kind.
template <typename T>
inline Enable_if<Real_type<T>(), T> comp_ellint_2(const T& k)
{
    // Algorithm:
    // ----------
    // Cody, W. J. Chebyshev approximations for the complete elliptic
    // integrals K and E. Mathematics of Computation, 1965, vol. 19,
    // pp. 105-112.
    //
    // Coefficients are taken from Table III.

    assert(k >= 0.0 && k < 1.0);
    const T eta = 1.0 - k * k;

    constexpr int n = 8;
    constexpr std::array<T, n> c = {
        4.4314718112155806e-1, 5.6805657874695358e-2, 2.1876220647186198e-2,
        1.2510592410844644e-2, 1.3034146073731432e-2, 1.5377102528552019e-2,
        7.3356174974290365e-3, 7.0980964089987229e-4};
    constexpr std::array<T, n> d = {
        2.4999999993617622e-1, 9.3749920249680113e-2, 5.8582839536559024e-2,
        4.2382807456947900e-2, 3.0302747728412848e-2, 1.5525129948040721e-2,
        3.4838679435896492e-3, 1.6427210797048025e-4};

    T res = 0.0;
    for (int i = 0; i < n; ++i) {
        res += d[i] * std::pow(eta, i + 1);
    }
    res *= std::log(1.0 / eta);
    for (int i = 0; i < n; ++i) {
        res += c[i] * std::pow(eta, i + 1);
    }
    return res + 1.0;
}

} // namespace Numlib

#endif // NUMLIB_MATH_SPECIAL_FUNCTIONS_H

