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

namespace num {

//------------------------------------------------------------------------------

// Equality comparable:

// Two matrices compare equal when they have the same elements. Comparison of
// matrices decribed by different slices is undefined behavior.

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
    Matrix<T, N> res = a;
    res += b;
    return res;
}

// Matrix subtraction:

template <typename T, std::size_t N>
Matrix<T, N> operator-(const Matrix<T, N>& a, const Matrix<T, N>& b)
{
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
