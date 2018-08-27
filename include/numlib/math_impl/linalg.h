// Copyright (c) 2018 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#ifndef NUMLIB_MATH_LINALG_H
#define NUMLIB_MATH_LINALG_H

#include <numlib/matrix.h>
#include <numlib/traits.h>
#include <algorithm>

namespace Numlib {

//------------------------------------------------------------------------------
//
// Find max, min, sum, and product of elements:

template <typename M>
inline Enable_if<Matrix_type<M>(), typename M::value_type> max(const M& vec)
{
    static_assert(M::order == 1, "bad rank for max(vec)");
    return *std::max_element(vec.begin(), vec.end());
}

template <typename M>
inline Enable_if<Matrix_type<M>(), Vec<typename M::value_type>>
max(const M& mat, std::size_t dim)
{
    static_assert(M::order == 2, "bad rank for max(mat)");
    assert(dim >= 0 && dim < M::order);

    Vec<Value_type<M>> result(mat.extent(dim));
    if (dim == 0) { // row
        for (std::size_t i = 0; i < mat.rows(); ++i) {
            result(i) = max(mat.row(i));
        }
    }
    else { // column
        for (std::size_t i = 0; i < mat.rows(); ++i) {
            result(i) = max(mat.column(i));
        }
    }
    return result;
}

} // namespace Numlib

#endif // NUMLIB_MATH_LINALG_H
