// Copyright (c) 2008-2010 Kent State University
// Copyright (c) 2011-2012 Texas A&M University
// Copyright (c) 2018 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#ifndef NUMLIB_BAND_MATRIX_OPERATIONS_H
#define NUMLIB_BAND_MATRIX_OPERATIONS_H

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

} // namespace Numlib

#endif // NUMLIB_BAND_MATRIX_OPERATIONS_H
