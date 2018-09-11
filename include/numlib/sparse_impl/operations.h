// Copyright (c) 2018 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#ifndef NUMLIB_SPARSE_OPERATIONS_H
#define NUMLIB_SPARSE_OPERATIONS_H

#include <numlib/matrix.h>
#include <vector>

namespace Numlib {

//------------------------------------------------------------------------------
//
// Convert formats:

// Gather a sparse full-storage vector into compressed form.
template <typename T>
Sparse_vector<T> gather(const Matrix<T, 1>& y)
{
    std::vector<T> val;
    std::vector<std::ptrdiff_t> loc;

    for (std::ptrdiff_t i = 0; i < y.size(); ++i) {
        if (y(i) != T{0}) {
            val.push_back(y(i));
            loc.push_back(i);
        }
    }
    return {val, loc};
}

// Scatter a sparse vector into full storage form.
template <typename T>
Matrix<T, 1> scatter(const Sparse_vector<T>& y)
{
    Matrix<T, 1> res(y.size());
    for (std::ptrdiff_t i = 0; i < res.size(); ++i) {
        res(i) = y(i);
    }
    return res;
}

} // namespace Numlib

#endif // NUMLIB_SPARSE_OPERATIONS_H
