// Copyright (c) 2008-2010 Kent State University
// Copyright (c) 2011-2012 Texas A&M University
// Copyright (c) 2018 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#ifndef NUMLIB_MATRIX_MATRIX_BASE_H
#define NUMLIB_MATRIX_MATRIX_BASE_H

#include <cassert>

namespace Numlib {

// Provides support for features common to both matrices and matrix reference.
template <typename T, std::size_t N>
class Matrix_base {
public:
    // Number of dimensions:
    static constexpr std::size_t order = N;

    using value_type = T;
    using size_type = std::ptrdiff_t;

    Matrix_base() = default;

    // Need a static_cast to avoid narrowing error:
    template <typename... Exts>
    explicit Matrix_base(Exts... exts) : desc{static_cast<size_type>(exts)...}
    {
    }

    explicit Matrix_base(const Matrix_slice<N>& s) : desc{s} {}

    // Move construction and assignment:
    Matrix_base(Matrix_base&&) = default;
    Matrix_base& operator=(Matrix_base&&) = default;

    // Copy construction and assignment:
    Matrix_base(const Matrix_base&) = default;
    Matrix_base& operator=(const Matrix_base&) = default;

    ~Matrix_base() = default;

    // Matrix rank.
    std::size_t rank() const { return N; }

    // Total number of elements:
    size_type size() const { return desc.size; }

    // Number of elements in the N'th dimension:
    size_type extent(size_type n) const
    {
        assert(n >= 0 && n < order);
        return desc.extents[n];
    }

    // Number of rows.
    size_type rows() const
    {
        assert(N >= 1);
        return desc.extents[0];
    }

    // Number of columns.
    size_type cols() const
    {
        assert(N >= 2);
        return desc.extents[1];
    }

    // The slice defining the shape of the matrix:
    const Matrix_slice<N>& shape() const { return desc; }
    const Matrix_slice<N>& descriptor() const { return desc; }

protected:
    Matrix_slice<N> desc;
};

} // namespace Numlib

#endif // NUMLIB_MATRIX_MATRIX_BASE_H
