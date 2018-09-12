// Copyright (c) 2018 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#ifndef NUMLIB_SPARSE_MATRIX_H
#define NUMLIB_SPARSE_MATRIX_H

#include <array>
#include <cassert>
#include <vector>
#include <numlib/traits.h>

namespace Numlib {

template <typename T>
class Sparse_matrix {
public:
    using value_type = T;
    using size_type = std::ptrdiff_t;
    using iterator = typename std::vector<T>::iterator;
    using const_iterator = typename std::vector<T>::const_iterator;

    Sparse_matrix() = default;

    // Copy semantics:
    Sparse_matrix(const Sparse_matrix&) = default;
    Sparse_matrix& operator=(const Sparse_matrix&) = default;

    // Move semantics:
    Sparse_matrix(Sparse_matrix&&) = default;
    Sparse_matrix& operator=(Sparse_matrix&&) = default;

    // Create from values, columns, and row indices.
    Sparse_matrix(size_type nr,
                  size_type nc,
                  const std::vector<T>& val,
                  const std::vector<size_type>& col,
                  const std::vector<size_type>& row)
        : elems(val), col_indx(col), row_ptr(row), extents{nr, nc}
    {
        assert(elems.size() == col_indx.size());
        assert(row_ptr.size() == narrow_cast<std::size_t>(nr + 1));
    }

    ~Sparse_matrix() = default;

private:
    std::vector<T> elems;
    std::vector<size_type> col_indx;
    std::vector<size_type> row_ptr;
    std::array<size_type, 2> extents;

    static const T zero;
};

// clang-format off
template <typename T>
const typename Sparse_matrix<T>::value_type 
Sparse_matrix<T>::zero = value_type{0};
// clang-format on

} // namespace Numlib

#endif // NUMLIB_SPARSE_MATRIX_H
