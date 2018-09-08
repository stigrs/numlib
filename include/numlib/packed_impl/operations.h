// Copyright (c) 2018 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#ifndef NUMLIB_PACKED_MATRIX_OPERATIONS_H
#define NUMLIB_PACKED_MATRIX_OPERATIONS_H

#include <algorithm>
#include <iomanip>
#include <iostream>

namespace Numlib {

//------------------------------------------------------------------------------
//
// Non-member functions returning properties of a packed matrix:

// Return size.
template <typename T, Uplo_scheme Uplo>
inline std::ptrdiff_t size(const Packed_matrix<T, Uplo>& ap)
{
    return ap.size();
}

// Return number of rows.
template <typename T, Uplo_scheme Uplo>
inline std::ptrdiff_t rows(const Packed_matrix<T, Uplo>& ap)
{
    return ap.rows();
}

// Return number of columns.
template <typename T, Uplo_scheme Uplo>
inline std::ptrdiff_t cols(const Packed_matrix<T, Uplo>& ap)
{
    return ap.cols();
}

// Return extent for a given dimension.
template <typename T, Uplo_scheme Uplo>
inline std::ptrdiff_t extent(const Packed_matrix<T, Uplo>& ap,
                             std::ptrdiff_t dim)
{
    assert(0 <= dim && dim < 2);
    return ap.extent(dim);
}

// Return UPLO scheme.
template <typename T, Uplo_scheme Uplo>
inline char uplo_scheme(const Packed_matrix<T, Uplo>& ap)
{
    return ap.uplo_scheme();
}

//------------------------------------------------------------------------------
//
// Equality comparable:

// Two matrices compare equal when they have the same elements. Comparison
// of matrices decribed by different slices is undefined behavior.

template <typename T, Uplo_scheme Uplo>
inline bool operator==(const Packed_matrix<T, Uplo>& a,
                       const Packed_matrix<T, Uplo>& b)
{
    assert(a.rows() == b.rows() && a.cols() == b.cols());
    return std::equal(a.begin(), a.end(), b.begin());
}

template <typename T, Uplo_scheme Uplo>
inline bool operator!=(const Packed_matrix<T, Uplo>& a,
                       const Packed_matrix<T, Uplo>& b)
{
    return !(a == b);
}

//------------------------------------------------------------------------------
//
// Binary arithmetic operations:

// Scalar addition:

template <typename T, Uplo_scheme Uplo>
inline Packed_matrix<T, Uplo> operator+(const Packed_matrix<T, Uplo>& a,
                                        const T& scalar)
{
    Packed_matrix<T, Uplo> res(a);
    return res += scalar;
}

template <typename T, Uplo_scheme Uplo>
inline Packed_matrix<T, Uplo> operator+(const T& scalar,
                                        const Packed_matrix<T, Uplo>& a)
{
    Packed_matrix<T, Uplo> res(a);
    return res += scalar;
}

// Scalar subtraction:

template <typename T, Uplo_scheme Uplo>
inline Packed_matrix<T, Uplo> operator-(const Packed_matrix<T, Uplo>& a,
                                        const T& scalar)
{
    Packed_matrix<T, Uplo> res(a);
    return res -= scalar;
}

// Scalar multiplication:

template <typename T, Uplo_scheme Uplo>
inline Packed_matrix<T, Uplo> operator*(const Packed_matrix<T, Uplo>& a,
                                        const T& scalar)
{
    Packed_matrix<T, Uplo> res(a);
    return res *= scalar;
}

template <typename T, Uplo_scheme Uplo>
inline Packed_matrix<T, Uplo> operator*(const T& scalar,
                                        const Packed_matrix<T, Uplo>& a)
{
    Packed_matrix<T, Uplo> res(a);
    return res *= scalar;
}

// Scalar division:

template <typename T, Uplo_scheme Uplo>
inline Packed_matrix<T, Uplo> operator/(const Packed_matrix<T, Uplo>& a,
                                        const T& scalar)
{
    Packed_matrix<T, Uplo> res(a);
    return res /= scalar;
}

// Scalar modulus:

template <typename T, Uplo_scheme Uplo>
inline Packed_matrix<T, Uplo> operator%(const Packed_matrix<T, Uplo>& a,
                                        const T& scalar)
{
    Packed_matrix<T, Uplo> res(a);
    return res %= scalar;
}

//------------------------------------------------------------------------------
//
// Output to stream:

template <typename T, Uplo_scheme Uplo>
std::ostream& operator<<(std::ostream& to, const Packed_matrix<T, Uplo>& ap)
{
    to << ap.rows() << " x " << ap.cols() << "\n[";
    for (std::ptrdiff_t i = 0; i < ap.rows(); ++i) {
        for (std::ptrdiff_t j = 0; j < ap.cols(); ++j) {
            to << std::setw(9) << ap(i, j) << " ";
        }
        if (i != ap.rows() - 1) {
            to << "\n ";
        }
    }
    to << "]\n";
    return to;
}

} // namespace Numlib

#endif // NUMLIB_PACKED_MATRIX_OPERATIONS_H
