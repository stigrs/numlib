// Copyright (c) 2018 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#ifndef NUMLIB_SPARSE_MATRIX_H
#define NUMLIB_SPARSE_MATRIX_H

#ifdef _OPENMP
#include <omp.h>
#endif

#include <algorithm>
#include <array>
#include <cassert>
#include <vector>
#include <numlib/traits.h>

namespace Numlib {

// Range-checked sparse matrix class.
//
// This class provides a basic framework for implementing sparse matrix
// methods that utilize the Intel Math Kernel Library.
//
// Note:
// - Elements are stored in the three array variation of the compressed
//   sparse row (CSR3) format.
// - Zero and one-based indexing are supported.
// - It is assumed that the sparse matrix is initialized with element indices
//   sorted in ascending order.
// - New elements are inserted so that the index order is preserved.
// - Size type is BLAS_INT to allow linking with Intel MKL.
//
template <typename T>
class Sparse_matrix {
public:
    using value_type = T;
    using size_type = BLAS_INT;
    using iterator = typename std::vector<T>::iterator;
    using const_iterator = typename std::vector<T>::const_iterator;

    Sparse_matrix() = default;

    // Copy semantics:
    Sparse_matrix(const Sparse_matrix&) = default;
    Sparse_matrix& operator=(const Sparse_matrix&) = default;

    // Move semantics:
    Sparse_matrix(Sparse_matrix&&) = default;
    Sparse_matrix& operator=(Sparse_matrix&&) = default;

    // Create from values, columns, and row indices:

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

    template <BLAS_INT n, BLAS_INT nnz>
    Sparse_matrix(size_type nr,
                  size_type nc,
                  const T (&val)[n],
                  const BLAS_INT (&col)[n],
                  const BLAS_INT (&row)[nnz]);

    ~Sparse_matrix() = default;

    // "Flat" element access:

    T* data() { return elems.data(); }
    const T* data() const { return elems.data(); }

    // Access underlying arrays:

    auto& values() { return elems; }
    const auto& values() const { return elems; }
    const auto& columns() const { return col_indx; }
    const auto& row_index() const { return row_ptr; }

    const auto& columns_zero_based() const { return col_indx; }
    const auto& row_index_zero_based() const { return row_ptr; }

    auto columns_one_based() const;
    auto row_index_one_based() const;

    // Properties:

    bool empty() const { return elems.empty(); }

    size_type size() const { return extents[0] * extents[1]; }
    size_type num_nonzero() const { return elems.size(); }
    size_type rows() const { return extents[0]; }
    size_type cols() const { return extents[1]; }
    size_type extent(size_type dim) const
    {
        assert(0 <= dim && dim < 2);
        return extents[dim];
    }

    // Subscripting.
    const T& operator()(size_type i, size_type j) const { return ref(i, j); }

    // Iterators:

    iterator begin() { return elems.begin(); }
    const_iterator begin() const { return elems.begin(); }

    iterator end() { return elems.end(); }
    const_iterator end() const { return elems.end(); }

    // Mutators:

    void swap(Sparse_matrix& m);
    void insert(size_type i, size_type j, const T& value);

    // Apply f(x) for every element x.
    template <typename F>
    Sparse_matrix& apply(F f);

    // Arithmetic operations:

    Sparse_matrix& operator*=(const T& value); // scalar multiplication
    Sparse_matrix& operator/=(const T& value); // scalar division

private:
    std::vector<T> elems;
    std::vector<size_type> col_indx;
    std::vector<size_type> row_ptr;
    std::array<size_type, 2> extents;

    static const T zero;

    const T& ref(size_type i, size_type j) const;
};

template <typename T>
template <BLAS_INT n, BLAS_INT nnz>
Sparse_matrix<T>::Sparse_matrix(size_type nr,
                                size_type nc,
                                const T (&val)[n],
                                const BLAS_INT (&col)[n],
                                const BLAS_INT (&row)[nnz])
    : elems(n), col_indx(n), row_ptr(nnz), extents{nr, nc}
{
    assert(row_ptr.size() == narrow_cast<std::size_t>(nr + 1));

    for (size_type i = 0; i < n; ++i) {
        elems[i] = val[i];
        col_indx[i] = col[i];
    }
    for (size_type i = 0; i < nnz; ++i) {
        row_ptr[i] = row[i];
    }
}

template <typename T>
auto Sparse_matrix<T>::columns_one_based() const
{
    auto result = col_indx;
    for (auto& i : result) {
        i += 1;
    }
    return result;
}

template <typename T>
auto Sparse_matrix<T>::row_index_one_based() const
{
    auto result = row_ptr;
    for (auto& i : result) {
        i += 1;
    }
    return result;
}

template <typename T>
void Sparse_matrix<T>::swap(Sparse_matrix& m)
{
    elems.swap(m.elems);
    col_indx.swap(m.col_indx);
    row_ptr.swap(m.row_ptr);
    std::swap(extents, m.extents);
}

template <typename T>
void Sparse_matrix<T>::insert(size_type i, size_type j, const T& value)
{
    if (ref(i, j) == zero) {
        auto pos = std::upper_bound(col_indx.begin() + row_ptr[i],
                                    col_indx.begin() + row_ptr[i + 1], j);
        size_type offset =
            narrow_cast<BLAS_INT>(std::distance(col_indx.begin(), pos));
        elems.insert(elems.begin() + offset, value);
        col_indx.insert(pos, j);
        for (std::size_t k = i + 1; k < row_ptr.size(); ++k) {
            row_ptr[k]++;
        }
    }
}

template <typename T>
template <typename F>
inline Sparse_matrix<T>& Sparse_matrix<T>::apply(F f)
{
#pragma omp parallel for
    for (std::size_t i = 0; i < elems.size(); ++i) {
        f(elems[i]);
    }
    return *this;
}

template <typename T>
inline Sparse_matrix<T>& Sparse_matrix<T>::operator*=(const T& value)
{
    return apply([&](T& a) { a *= value; });
}

template <typename T>
inline Sparse_matrix<T>& Sparse_matrix<T>::operator/=(const T& value)
{
    return apply([&](T& a) { a /= value; });
}

template <typename T>
const T& Sparse_matrix<T>::ref(size_type i, size_type j) const
{
    for (size_type k = row_ptr[i]; k < row_ptr[i + 1]; ++k) {
        if (col_indx[k] == j) {
            return elems[k];
        }
    }
    return zero;
}

// clang-format off
template <typename T>
const typename Sparse_matrix<T>::value_type 
Sparse_matrix<T>::zero = value_type{0};
// clang-format on

} // namespace Numlib

#endif // NUMLIB_SPARSE_MATRIX_H
