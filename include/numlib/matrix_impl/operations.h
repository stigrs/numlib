// Copyright (c) 2008-2010 Kent State University
// Copyright (c) 2011-2012 Texas A&M University
// Copyright (c) Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#ifndef NUMLIB_MATRIX_OPERATIONS_H
#define NUMLIB_MATRIX_OPERATIONS_H

namespace Numlib {

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

} // namespace Numlib

#endif // NUMLIB_MATRIX_OPERATIONS_H
