// Copyright (c) 2018 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#ifndef NUMLIB_MATH_CALCULUS_H
#define NUMLIB_MATH_CALCULUS_H

#include <numlib/matrix.h>
#include <numlib/traits.h>
#include <array>
#include <cassert>
#include <cmath>
#include <functional>
#include <limits>
#include <string>

#ifdef ENABLE_QUADPACK
#include <numlib/math_impl/quadpack.h>
#endif

namespace Numlib {

//------------------------------------------------------------------------------
//
// Numerical derivation:

// Compute the numerical first derivative of the function f(x).
inline double dfdx(std::function<double(double)> f, double x)
{
    auto eps = std::numeric_limits<double>::epsilon();
    auto h = std::pow(eps, 1.0 / 3.0) * x; // Numerical recipes
    return (f(x + h) - f(x - h)) / (2.0 * h);
}

//------------------------------------------------------------------------------
//
// Numerical integration:

// Integrate function values over a non-uniform grid using the
// Trapezoidal rule.
double trapz(double xlo, double xup, const Vec<double>& y);

// Return tabulated roots and weights for a Gauss-Legendre quadrature
// of order n.
template <int N = 8>
void gauss_legendre(std::array<double, N>& roots,
                    std::array<double, N>& weights,
                    double a = -1.0,
                    double b = 1.0)
{
    static_assert(N == 5 || N == 8 || N == 16,
                  "bad order for Gauss-Legendre quadrature");

    // 5-point:
    std::array<double, 5> x5{
        0.00000000000000000, -0.5384693101056831, 0.5384693101056831,
        -0.9061798459386640, 0.9061798459386640,
    };
    std::array<double, 5> w5{
        0.5688888888888889, 0.4786286704993665, 0.4786286704993665,
        0.2369268850561891, 0.2369268850561891,
    };
    std::array<double, 8> x8{-0.1834346424956498, 0.1834346424956498,
                             -0.5255324099163290, 0.5255324099163290,
                             -0.7966664774136267, 0.7966664774136267,
                             -0.9602898564975363, 0.9602898564975363};
    std::array<double, 8> w8{0.3626837833783620, 0.3626837833783620,
                             0.3137066458778873, 0.3137066458778873,
                             0.2223810344533745, 0.2223810344533745,
                             0.1012285362903763, 0.1012285362903763};

    std::array<double, 16> x16{
        -0.0950125098376374, 0.0950125098376374,  -0.2816035507792589,
        0.2816035507792589,  -0.4580167776572274, 0.4580167776572274,
        -0.6178762444026438, 0.6178762444026438,  -0.7554044083550030,
        0.7554044083550030,  -0.8656312023878318, 0.8656312023878318,
        -0.9445750230732326, 0.9445750230732326,  -0.9894009349916499,
        0.9894009349916499};
    std::array<double, 16> w16{
        0.1894506104550685, 0.1894506104550685, 0.1826034150449236,
        0.1826034150449236, 0.1691565193950025, 0.1691565193950025,
        0.1495959888165767, 0.1495959888165767, 0.1246289712555339,
        0.1246289712555339, 0.0951585116824928, 0.0951585116824928,
        0.0622535239386479, 0.0622535239386479, 0.0271524594117541,
        0.0271524594117541};

    switch (N) {
    case 16:
        std::copy(x16.begin(), x16.end(), roots.begin());
        std::copy(w16.begin(), w16.end(), weights.begin());
        break;
    case 5:
        std::copy(x5.begin(), x5.end(), roots.begin());
        std::copy(w5.begin(), w5.end(), weights.begin());
        break;
    case 8:
    default:
        std::copy(x8.begin(), x8.end(), roots.begin());
        std::copy(w8.begin(), w8.end(), weights.begin());
    }

    // Change of interval:
    for (int i = 0; i < N; ++i) {
        roots[i] = 0.5 * (b - a) * roots[i] + 0.5 * (a + b);
        weights[i] *= 0.5 * (b - a);
    }
}

// Integrate function from a to b using a Gauss-Legendre quadrature
// of order N.
template <int N = 8>
double quad(std::function<double(double)> f, double a, double b)
{
    static_assert(N == 5 || N == 8 || N == 16,
                  "bad order for Gauss-Legendre quadrature");

    std::array<double, N> x;
    std::array<double, N> w;

    gauss_legendre<N>(x, w, a, b);

    double res = 0.0;
    for (int i = 0; i < N; ++i) {
        res += w[i] * f(x[i]);
    }
    return res;
}

#ifdef ENABLE_QUADPACK
// Wrapper to the QAGS subroutine from QUADPACK.
inline double qags(quadpack_fptr f,
                   double a,
                   double b,
                   double epsabs = 1.0e-15,
                   double epsrel = 1.0e-12,
                   int limit = 1000)
{
    int neval;
    int ier;
    int lenw = limit * 4;
    int last;

    Vec<int> iwork(limit);
    Vec<double> work(lenw);

    double result;
    double abserr;

    dqags_(f, a, b, epsabs, epsrel, result, abserr, neval, ier, limit, lenw,
           last, iwork.data(), work.data());
    if (ier != 0) {
        throw Math_error("dqags_ failed with error " + std::to_string(ier));
    }
    return result;
}

// Wrapper to the QAGI subroutine from QUADPACK.
inline double qagi(quadpack_fptr f,
                   double bound,
                   int inf,
                   double epsabs = 1.0e-15,
                   double epsrel = 1.0e-12,
                   int limit = 1000)
{
    assert(inf == 1 || inf == -1 || inf == 2); // see description for dqagi

    int neval;
    int ier;
    int lenw = limit * 4;
    int last;

    Vec<int> iwork(limit);
    Vec<double> work(lenw);

    double result;
    double abserr;

    dqagi_(f, bound, inf, epsabs, epsrel, result, abserr, neval, ier, limit,
           lenw, last, iwork.data(), work.data());
    if (ier != 0) {
        throw Math_error("dqagi_ failed with error " + std::to_string(ier));
    }
    return result;
}
#endif // ENABLE_QUADPACK

//------------------------------------------------------------------------------
//
// Solve ordinary differential equation (ODE):

// Fourth-order Runge-Kutta method.
inline double
rk4(std::function<double(double, double)> f, double dx, double x, double y)
{
    double k1 = dx * f(x, y);
    double k2 = dx * f(x + dx / 2.0, y + k1 / 2.0);
    double k3 = dx * f(x + dx / 2.0, y + k2 / 2.0);
    double k4 = dx * f(x + dx, y + k3);

    return y + (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0;
}

// Fourth-order Runge-Kutta method.
void rk4(std::function<Vec<double>(double t, const Vec<double>&)> f,
         Vec<double>& y,
         double t0,
         double t1,
         double dt);

} // namespace Numlib

#endif // NUMLIB_MATH_CALCULUS_H
