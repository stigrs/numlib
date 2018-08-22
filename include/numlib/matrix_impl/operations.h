// Copyright (c) 2008-2010 Kent State University
// Copyright (c) 2011-2012 Texas A&M University
// Copyright (c) 2018 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#ifndef NUMLIB_MATRIX_OPERATIONS_H
#define NUMLIB_MATRIX_OPERATIONS_H

#include <algorithm>
#include <random>

namespace num {

//------------------------------------------------------------------------------
//
// The following operations are defined for all Matrix types:

// Return matrix rank.
template <typename M>
inline Enable_if<Matrix_type<M>(), std::size_t> rank(const M& m)
{
    return m.rank();
}

// Return matrix size.
template <typename M>
inline Enable_if<Matrix_type<M>(), std::size_t> size(const M& m)
{
    return m.size();
}

// Return number of rows.
template <typename M>
inline Enable_if<Matrix_type<M>(), std::size_t> rows(const M& m)
{
    static_assert(0 < M::order, "");
    return m.extent(0);
}

// Return number of columns.
template <typename M>
inline Enable_if<Matrix_type<M>(), std::size_t> cols(const M& m)
{
    static_assert(0 < M::order, "");
    return m.extent(1);
}

// Create matrix of zeros.
template <typename M, typename... Args>
inline Enable_if<Matrix_type<M>(), M> zeros(Args... args)
{
    using value_type = typename M::value_type;

    assert(M::order == sizeof...(args));
    M res(args...);
    res = value_type{0};
    return res;
}

// Create matrix of ones.
template <typename M, typename... Args>
inline Enable_if<Matrix_type<M>(), M> ones(Args... args)
{
    using value_type = typename M::value_type;

    assert(M::order == sizeof...(args));
    M res(args...);
    res = value_type{1};
    return res;
}

// Create a random matrix from a normal distribution with zero mean and unit
// variance.
template <typename M, typename... Args>
inline Enable_if<Matrix_type<M>() && Real_type<typename M::value_type>(), M>
randn(Args... args)
{
    assert(M::order == sizeof...(args));
    M res(args...);

    std::random_device rd{};
    std::mt19937_64 gen{rd()};
    std::normal_distribution<> nd{};

    for (auto& x : res) {
        x = nd(gen);
    }
    return res;
}

// Create a random matrix from a uniform real distribution on the
// interval [0, 1).
template <typename M, typename... Args>
inline Enable_if<Matrix_type<M>() && Real_type<typename M::value_type>(), M>
randu(Args... args)
{
    assert(M::order == sizeof...(args));
    M res(args...);

    std::random_device rd{};
    std::mt19937_64 gen{rd()};

    std::uniform_real_distribution<> ur{};
    for (auto& x : res) {
        x = ur(gen);
    }
    return res;
}

// Create a random matrix from a uniform integer distribution on the
// interval [0, 1].
template <typename M, typename... Args>
inline Enable_if<Matrix_type<M>() && Integer_type<typename M::value_type>(), M>
randi(Args... args)
{
    assert(M::order == sizeof...(args));
    M res(args...);

    std::random_device rd{};
    std::mt19937_64 gen{rd()};

    std::uniform_int_distribution<> ui{};
    for (auto& x : res) {
        x = ui(gen);
    }
    return res;
}

//------------------------------------------------------------------------------
//
// Equality comparable:

// Two matrices compare equal when they have the same elements. Comparison
// of matrices decribed by different slices is undefined behavior.

template <typename M1, typename M2>
inline Enable_if<Matrix_type<M1>() && Matrix_type<M2>(), bool>
operator==(const M1& a, const M2& b)
{
    assert(same_extents(a, b));
    return std::equal(a.begin(), a.end(), b.begin());
}

template <typename M1, typename M2>
inline Enable_if<Matrix_type<M1>() && Matrix_type<M2>(), bool>
operator!=(const M1& a, const M2& b)
{
    return !(a == b);
}

//------------------------------------------------------------------------------

// Binary arithmetic operations:

// Matrix addition:

template <typename T, std::size_t N>
Matrix<T, N> operator+(const Matrix<T, N>& a, const Matrix<T, N>& b)
{
    assert(same_extents(a, b));

    Matrix<T, N> res = a;
    res += b;
    return res;
}

// Matrix subtraction:

template <typename T, std::size_t N>
Matrix<T, N> operator-(const Matrix<T, N>& a, const Matrix<T, N>& b)
{
    assert(same_extents(a, b));

    Matrix<T, N> res = a;
    res -= b;
    return res;
}

// Scalar multiplication:

template <typename T, std::size_t N>
Matrix<T, N> operator*(const Matrix<T, N>& a, const T& scalar)
{
    Matrix<T, N> res = a;
    res *= scalar;
    return res;
}

template <typename T, std::size_t N>
Matrix<T, N> operator*(const T& scalar, const Matrix<T, N>& a)
{
    Matrix<T, N> res = a;
    res *= scalar;
    return res;
}

} // namespace num

#endif // NUMLIB_MATRIX_OPERATIONS_H
