// Copyright (c) 2018 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#ifndef NUMLIB_BAND_MATRIX_OPERATIONS_H
#define NUMLIB_BAND_MATRIX_OPERATIONS_H

#include <iomanip>
#include <iostream>
#include <algorithm>

namespace Numlib {

//------------------------------------------------------------------------------

// Non-member functions returning properties of a band matrix:

// Return band matrix size.
template <typename T>
inline std::ptrdiff_t size(const Band_matrix<T>& ab)
{
    return ab.size();
}

// Return number of rows.
template <typename T>
inline std::ptrdiff_t rows(const Band_matrix<T>& ab)
{
    return ab.rows();
}

// Return number of columns.
template <typename T>
inline std::ptrdiff_t cols(const Band_matrix<T>& ab)
{
    return ab.cols();
}

// Return extent for a given dimension.
template <typename T>
inline std::ptrdiff_t extent(const Band_matrix<T>& ab, std::ptrdiff_t dim)
{
    assert(0 <= dim && dim < 2);
    return ab.extent(dim);
}

// Return bandwidth.
template <typename T>
inline std::ptrdiff_t bwidth(const Band_matrix<T>& ab, std::ptrdiff_t uplo)
{
    assert(0 <= uplo && uplo < 2);
    if (uplo == 0) {
        return ab.lower();
    }
    else {
        return ab.upper();
    }
}

//------------------------------------------------------------------------------
//
// Equality comparable:

// Two band matrices compare equal when they have the same elements.

template <typename T>
inline bool operator==(const Band_matrix<T>& a, const Band_matrix<T>& b)
{
    assert(a.rows() == b.rows() && a.cols() == b.cols());
    return std::equal(a.begin(), a.end(), b.begin());
}

template <typename T>
inline bool operator!=(const Band_matrix<T>& a, const Band_matrix<T>& b)
{
    return !(a == b);
}

//------------------------------------------------------------------------------
//
// Binary arithmetic operations:

// Scalar addition:

template <typename T>
inline Band_matrix<T> operator+(const Band_matrix<T>& a, const T& scalar)
{
    Band_matrix<T> res(a);
    return res += scalar;
}

template <typename T>
inline Band_matrix<T> operator+(const T& scalar, const Band_matrix<T>& a)
{
    Band_matrix<T> res(a);
    return res += scalar;
}

// Scalar subtraction:

template <typename T>
inline Band_matrix<T> operator-(const Band_matrix<T>& a, const T& scalar)
{
    Band_matrix<T> res(a);
    return res -= scalar;
}

// Scalar multiplication:

template <typename T>
inline Band_matrix<T> operator*(const Band_matrix<T>& a, const T& scalar)
{
    Band_matrix<T> res(a);
    return res *= scalar;
}

template <typename T>
inline Band_matrix<T> operator*(const T& scalar, const Band_matrix<T>& a)
{
    Band_matrix<T> res(a);
    return res *= scalar;
}

// Scalar division:

template <typename T>
inline Band_matrix<T> operator/(const Band_matrix<T>& a, const T& scalar)
{
    Band_matrix<T> res(a);
    return res /= scalar;
}

// Scalar modulus:

template <typename T>
inline Band_matrix<T> operator%(const Band_matrix<T>& a, const T& scalar)
{
    Band_matrix<T> res(a);
    return res %= scalar;
}

//------------------------------------------------------------------------------
//
// Output to stream:

template <typename T>
std::ostream& operator<<(std::ostream& to, const Band_matrix<T>& ab)
{
    to << ab.rows() << " x " << ab.cols() << "\n[";
    for (std::ptrdiff_t i = 0; i < ab.rows(); ++i) {
        for (std::ptrdiff_t j = 0; j < ab.cols(); ++j) {
            to << std::setw(9) << ab(i, j) << " ";
        }
        if (i != ab.rows() - 1) {
            to << "\n ";
        }
    }
    to << "]\n";
    return to;
}

} // namespace Numlib

#endif // NUMLIB_BAND_MATRIX_OPERATIONS_H
