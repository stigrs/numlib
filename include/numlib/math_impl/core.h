// Copyright (c) 2018 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#ifndef NUMLIB_MATH_CORE_H
#define NUMLIB_MATH_CORE_H

#include <numlib/traits.h>
#include <numlib/matrix.h>
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
    return static_cast<T>(std::round(x));
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

//------------------------------------------------------------------------------
//
// Miscellaneous element-wise functions:

template <typename T, std::size_t N>
inline Matrix<T, N> abs(const Matrix<T, N>& m)
{
    Matrix<T, N> res(m);
    res.apply([](T& x) { x = std::abs(x); });
    return res;
}

template <typename T, std::size_t N>
inline Enable_if<Real_type<T>(), Matrix<T, N>> pow(const Matrix<T, N>& m,
                                                   const T& val)
{
    Matrix<T, N> res(m);
    res.apply([](T& x, const T& p) { x = std::pow(x, p); }, val);
    return res;
}

template <typename T, std::size_t N>
inline Enable_if<Real_type<T>(), Matrix<T, N>> sqrt(const Matrix<T, N>& m)
{
    Matrix<T, N> res(m);
    res.apply([](T& x) { x = std::sqrt(x); });
    return res;
}

template <typename T, std::size_t N>
inline Enable_if<Real_type<T>(), Matrix<T, N>> cbrt(const Matrix<T, N>& m)
{
    Matrix<T, N> res(m);
    res.apply([](T& x) { x = std::cbrt(x); });
    return res;
}

template <typename T, std::size_t N>
inline Enable_if<Real_type<T>(), Matrix<T, N>> exp(const Matrix<T, N>& m)
{
    Matrix<T, N> res(m);
    res.apply([](T& x) { x = std::exp(x); });
    return res;
}

template <typename T, std::size_t N>
inline Enable_if<Real_type<T>(), Matrix<T, N>> log(const Matrix<T, N>& m)
{
    Matrix<T, N> res(m);
    res.apply([](T& x) { x = std::log(x); });
    return res;
}

template <typename T, std::size_t N>
inline Enable_if<Real_type<T>(), Matrix<T, N>> erf(const Matrix<T, N>& m)
{
    Matrix<T, N> res(m);
    res.apply([](T& x) { x = std::erf(x); });
    return res;
}

template <typename T, std::size_t N>
inline Enable_if<Real_type<T>(), Matrix<T, N>> erfc(const Matrix<T, N>& m)
{
    Matrix<T, N> res(m);
    res.apply([](T& x) { x = std::erfc(x); });
    return res;
}

template <typename T, std::size_t N>
inline Enable_if<Real_type<T>(), Matrix<T, N>> tgamma(const Matrix<T, N>& m)
{
    Matrix<T, N> res(m);
    res.apply([](T& x) { x = std::tgamma(x); });
    return res;
}

template <typename T, std::size_t N>
inline Enable_if<Real_type<T>(), Matrix<T, N>> lgamma(const Matrix<T, N>& m)
{
    Matrix<T, N> res(m);
    res.apply([](T& x) { x = std::lgamma(x); });
    return res;
}

template <typename T, std::size_t N>
inline Enable_if<Real_type<T>(), Matrix<T, N>> sin(const Matrix<T, N>& m)
{
    Matrix<T, N> res(m);
    res.apply([](T& x) { x = std::sin(x); });
    return res;
}

template <typename T, std::size_t N>
inline Enable_if<Real_type<T>(), Matrix<T, N>> cos(const Matrix<T, N>& m)
{
    Matrix<T, N> res(m);
    res.apply([](T& x) { x = std::cos(x); });
    return res;
}

template <typename T, std::size_t N>
inline Enable_if<Real_type<T>(), Matrix<T, N>> tan(const Matrix<T, N>& m)
{
    Matrix<T, N> res(m);
    res.apply([](T& x) { x = std::tan(x); });
    return res;
}

template <typename T, std::size_t N>
inline Enable_if<Real_type<T>(), Matrix<T, N>> asin(const Matrix<T, N>& m)
{
    Matrix<T, N> res(m);
    res.apply([](T& x) { x = std::asin(x); });
    return res;
}

template <typename T, std::size_t N>
inline Enable_if<Real_type<T>(), Matrix<T, N>> acos(const Matrix<T, N>& m)
{
    Matrix<T, N> res(m);
    res.apply([](T& x) { x = std::acos(x); });
    return res;
}

template <typename T, std::size_t N>
inline Enable_if<Real_type<T>(), Matrix<T, N>> atan(const Matrix<T, N>& m)
{
    Matrix<T, N> res(m);
    res.apply([](T& x) { x = std::atan(x); });
    return res;
}

template <typename T, std::size_t N>
inline Enable_if<Real_type<T>(), Matrix<T, N>> sinh(const Matrix<T, N>& m)
{
    Matrix<T, N> res(m);
    res.apply([](T& x) { x = std::sinh(x); });
    return res;
}

template <typename T, std::size_t N>
inline Enable_if<Real_type<T>(), Matrix<T, N>> cosh(const Matrix<T, N>& m)
{
    Matrix<T, N> res(m);
    res.apply([](T& x) { x = std::cosh(x); });
    return res;
}

template <typename T, std::size_t N>
inline Enable_if<Real_type<T>(), Matrix<T, N>> tanh(const Matrix<T, N>& m)
{
    Matrix<T, N> res(m);
    res.apply([](T& x) { x = std::tanh(x); });
    return res;
}

template <typename T, std::size_t N>
inline Enable_if<Real_type<T>(), Matrix<T, N>> asinh(const Matrix<T, N>& m)
{
    Matrix<T, N> res(m);
    res.apply([](T& x) { x = std::asinh(x); });
    return res;
}

template <typename T, std::size_t N>
inline Enable_if<Real_type<T>(), Matrix<T, N>> acosh(const Matrix<T, N>& m)
{
    Matrix<T, N> res(m);
    res.apply([](T& x) { x = std::acosh(x); });
    return res;
}

template <typename T, std::size_t N>
inline Enable_if<Real_type<T>(), Matrix<T, N>> atanh(const Matrix<T, N>& m)
{
    Matrix<T, N> res(m);
    res.apply([](T& x) { x = std::atanh(x); });
    return res;
}

} // namespace Numlib

#endif // NUMLIB_MATH_CORE_H
