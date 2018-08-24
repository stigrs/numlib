// Copyright (c) 2018 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#ifndef NUMLIB_MATH_CALCULUS_H
#define NUMLIB_MATH_CALCULUS_H

#include <numlib/matrix.h>
#include <numlib/traits.h>
#include <cmath>
#include <functional>
#include <limits>

namespace Numlib {

namespace Math {

    //--------------------------------------------------------------------------
    //
    // Numerical derivation:

    // Compute the numerical first derivative of the function f(x).
    inline double dfdx(std::function<double(double)> f, double x)
    {
        auto eps = std::numeric_limits<double>::epsilon();
        auto h = std::pow(eps, 1.0 / 3.0) * x; // Numerical recipes
        return (f(x + h) - f(x - h)) / (2.0 * h);
    }

    //--------------------------------------------------------------------------
    //
    // Numerical integration:

    // Integrate function values over a non-uniform grid using the
    // Trapezoidal rule.
    double trapezoidal(double xlo, double xup, const Vec<double>& y);

} // namespace Math

} // namespace Numlib

#endif // NUMLIB_MATH_CALCULUS_H
