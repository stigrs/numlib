// Copyright (c) 2018 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#ifndef NUMLIB_MATH_CALCULUS_H
#define NUMLIB_MATH_CALCULUS_H

#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable : 4996) // caused by boost/numeric/odeint.hpp
#pragma warning(disable : 4127) // caused by boost/numeric/odeint.hpp
#endif

#ifdef __APPLE__
#pragma clang system_header
#endif
#include <boost/numeric/odeint.hpp>

#ifdef _MSC_VER
#pragma warning(pop)
#endif

#include <numlib/matrix.h>
#include <numlib/traits.h>
#include <array>
#include <cassert>
#include <cmath>
#include <functional>
#include <limits>
#include <string>

#include <numlib/math_impl/cquadpack.h>

namespace boost {
namespace numeric {
    namespace odeint {

        template <typename T>
        struct is_resizeable<Numlib::Vec<T>> {
            using type = boost::true_type;
            static const bool value = type::value;
        };
    } // namespace odeint
} // namespace numeric
} // namespace boost

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

// Wrapper to the QAGS subroutine from QUADPACK.
inline double qags(dq_function_type f,
                   double a,
                   double b,
                   double epsabs = 1.0e-15,
                   double epsrel = 1.0e-12)
{
    int neval;
    int ier;

    double res;
    double abserr;

    res = cquadpack_dqags(f, a, b, epsabs, epsrel, &abserr, &neval, &ier,
                          nullptr);
    if (ier != 0) {
        throw Math_error("dqags failed with error " + std::to_string(ier));
    }
    return res;
}

// Wrapper to the QAGI subroutine from QUADPACK.
inline double qagi(dq_function_type f,
                   double bound,
                   int inf,
                   double epsabs = 1.0e-15,
                   double epsrel = 1.0e-12)
{
    assert(inf == 1 || inf == -1 || inf == 2); // see description for dqagi

    int neval;
    int ier;

    double res;
    double abserr;

    res = cquadpack_dqagi(f, bound, inf, epsabs, epsrel, &abserr, &neval, &ier,
                          nullptr);
    if (ier != 0) {
        throw Math_error("dqagi failed with error " + std::to_string(ier));
    }
    return res;
}

//------------------------------------------------------------------------------
//
// Ordinary Differential Equation (ODE) Solvers:

// Solve an initial value problem for a system of ODEs.
//
// This is a wrapper to the Runge-Kutta DOPRI5 method provided by Boost Odeint.
//
inline void
solve_ivp(std::function<void(const Vec<double>&, Vec<double>&, const double)> f,
          Vec<double>& y,
          double& t0,
          double t1,
          double dt = 0.0,
          double atol = 1.0e-6,
          double rtol = 1.0e-6)
{
    using namespace boost::numeric::odeint;
    using state_t = Vec<double>;
    using error_stepper_t = runge_kutta_dopri5<state_t>;

    BOOST_STATIC_ASSERT(is_resizeable<state_t>::value == true);

    if (dt == 0.0) {
        dt = std::max(0.01 * (t1 - t0),
                      100.0 * std::numeric_limits<double>::epsilon());
    }
    integrate_adaptive(make_controlled<error_stepper_t>(atol, rtol), f, y, t0,
                       t1, dt);
    t0 = t1;
}

} // namespace Numlib

#endif // NUMLIB_MATH_CALCULUS_H
