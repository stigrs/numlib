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

    constexpr int n = 10;
    constexpr std::array<T, n> a = {
        9.6573590280856255384e-2, 3.0885146271305189866e-2,
        1.4938013532687165242e-2, 8.7898018745550646778e-3,
        6.1796274460533176084e-3, 6.8479092826245051197e-3,
        9.8489293221768937682e-3, 8.0030039806499853708e-3,
        2.2966348983969586869e-3, 1.3930878570066467279e-4};
    constexpr std::array<T, n> b = {
        1.2499999999990808051e-1, 7.0312499739038352054e-2,
        4.8828041906862397978e-2, 3.7377739758623604144e-2,
        3.0124849012898930266e-2, 2.3931913323110790077e-2,
        1.5530941631977203877e-2, 5.9739042991554291551e-3,
        9.2155463496324984638e-4, 2.9700280966555612066e-5};

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

    constexpr int n = 10;
    constexpr std::array<T, n> c = {
        4.4314718056088952648e-1, 5.6805194567559156648e-2,
        2.1831811676130481568e-2, 1.1569595745295402175e-2,
        7.5950934225594322802e-3, 7.8204040609595541727e-3,
        1.0770635039866455473e-2, 8.6384421736040744302e-3,
        2.4685033304607227339e-3, 1.4946621757181326771e-4};
    constexpr std::array<T, n> d = {
        2.4999999999990177208e-1, 9.3749999721203140658e-2,
        5.8593661255531491732e-2, 4.2717890547383095644e-2,
        3.3478943665761626232e-2, 2.6145014700313878932e-2,
        1.6804023346363384981e-2, 6.4321465864383017666e-3,
        9.8983328462253847867e-4, 3.1859195655501571800e-5};

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

