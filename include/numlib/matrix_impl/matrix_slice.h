// Copyright (c) 2008-2010 Kent State University
// Copyright (c) 2011-2012 Texas A&M University
// Copyright (c) 2018 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#ifndef NUMLIB_DENSE_MATRIX_SLICE_H
#define NUMLIB_DENSE_MATRIX_SLICE_H

#include <algorithm>
#include <array>
#include <cassert>
#include <initializer_list>
#include <iterator>

namespace Numlib {

// Matrix slice.
//
// A matrix slice specifies the N-dimensional matrix properties of a
// contigous region of memory. The slice is described by three parameters:
//
// - Starting offset
// - Sequence of extents
// - Sequence of strides
//
template <std::size_t N>
struct Matrix_slice {
    // Empty matrix:
    Matrix_slice() : size{0}, start{0} {}

    // Copy semantics:
    Matrix_slice(const Matrix_slice&) = default;
    Matrix_slice& operator=(const Matrix_slice&) = default;

    // Starting offset and extents:
    Matrix_slice(Index s, std::initializer_list<Index> exts);

    // Starting offset, extents, and strides:
    Matrix_slice(Index s,
                 std::initializer_list<Index> exts,
                 std::initializer_list<Index> strs);

    // N extents:
    template <typename... Dims>
    Matrix_slice(Dims... dims);

    // Calculate index from a set of subscripts:
    template <typename... Dims>
#ifdef _MSC_VER // Workaround for internal compiler error in VS 2017
    Index operator()(Dims... dims) const;
#else
    Enable_if<All(Convertible<Dims, Index>()...), Index>
    operator()(Dims... dims) const;
#endif // _MSC_VER

    // Calculate offset given a range.
    template <typename R>
    Index offset(R&& range) const;

    Index size;                   // total number of elements
    Index start;                  // starting offset
    std::array<Index, N> extents; // number of elements in each dim
    std::array<Index, N> strides; // offsets between elements
};

template <std::size_t N>
Matrix_slice<N>::Matrix_slice(Index s, std::initializer_list<Index> exts)
    : start(s)
{
    assert(exts.size() == N);
    std::copy(exts.begin(), exts.end(), extents.begin());
    Matrix_impl::compute_strides(*this);
}

template <std::size_t N>
Matrix_slice<N>::Matrix_slice(Index s,
                              std::initializer_list<Index> exts,
                              std::initializer_list<Index> strs)
    : start(s)
{
    assert(exts.size() == N);
    std::copy(exts.begin(), exts.end(), extents.begin());
    std::copy(strs.begin(), strs.end(), strides.begin());
    size = Matrix_impl::compute_size(extents);
}

template <std::size_t N>
template <typename... Dims>
Matrix_slice<N>::Matrix_slice(Dims... dims) : start{0}
{
    static_assert(sizeof...(Dims) == N,
                  "Matrix_slice<N>::Matrix_slice(Dims...): dimension mismatch");
    Index args[N]{Index(dims)...};
    std::copy(std::begin(args), std::end(args), extents.begin());
    Matrix_impl::compute_strides(*this);
}

template <std::size_t N>
template <typename... Dims>
#ifdef _MSC_VER // Workaround for internal compiler error in VS 2017
inline Index
#else
inline Enable_if<All(Convertible<Dims, Index>()...), Index>
#endif // _MSC_VER
Matrix_slice<N>::operator()(Dims... dims) const
{
    static_assert(sizeof...(Dims) == N,
                  "Matrix_slice<N>::operator(): dimension mismatch");
    Index args[N]{Index(dims)...};
    return start +
           std::inner_product(args, args + N, strides.begin(), Index{0});
}

template <std::size_t N>
template <typename R>
inline Index Matrix_slice<N>::offset(R&& range) const
{
    constexpr Index zero = 0;
    return start + std::inner_product(strides.begin(), strides.end(),
                                      std::begin(range), zero);
}

//------------------------------------------------------------------------------

// Specializations to allow optimizations:

// Matrix_slice to describe one-dimensional matrix (vector).
template <>
struct Matrix_slice<1> {
    // Empty matrix:
    Matrix_slice() = default;

    // Copy semantics:
    Matrix_slice(const Matrix_slice&) = default;
    Matrix_slice& operator=(const Matrix_slice&) = default;

    // Starting offset and extents:
    Matrix_slice(Index s, Index exts)
    {
        start = s;
        extents[0] = exts;
        strides[0] = 1;
        size = exts;
    }

    // Starting offset, extents, and strides:
    Matrix_slice(Index s, Index exts, Index strs)
    {
        start = s;
        extents[0] = exts;
        strides[0] = strs;
        size = exts * strs;
    }

    // N extents:
    Matrix_slice(Index exts)
    {
        start = 0;
        extents[0] = exts;
        strides[0] = 1;
        size = exts;
    }

    // Calculate index from a set of subscripts:
    Index operator()(Index i) const { return start + i * strides[0]; }

    // Calculate offset given a range.
    template <typename R>
    Index offset(R&& range) const
    {
        constexpr Index zero = 0;
        return start + std::inner_product(strides.begin(), strides.end(),
                                          std::begin(range), zero);
    }

    Index size;                   // total number of elements
    Index start;                  // starting offset
    std::array<Index, 1> extents; // number of elements in each dimension
    std::array<Index, 1> strides; // offsets between elements in each dim
};

// Matrix_slice to describe two-dimensional matrix (matrix).
template <>
struct Matrix_slice<2> {
    // Empty matrix:
    Matrix_slice() = default;

    // Copy semantics:
    Matrix_slice(const Matrix_slice&) = default;
    Matrix_slice& operator=(const Matrix_slice&) = default;

    // Starting offset and extents:
    Matrix_slice(Index s, std::initializer_list<Index> exts) : start(s)
    {
        assert(exts.size() == 2);
        std::copy(exts.begin(), exts.end(), extents.begin());
        Matrix_impl::compute_strides(*this);
    }

    // Starting offset, extents, and strides:
    Matrix_slice(Index s,
                 std::initializer_list<Index> exts,
                 std::initializer_list<Index> strs)
        : start(s)
    {
        assert(exts.size() == 2);
        std::copy(exts.begin(), exts.end(), extents.begin());
        std::copy(strs.begin(), strs.end(), strides.begin());
        size = Matrix_impl::compute_size(extents);
    }

    // N extents:
    Matrix_slice(Index nr, Index nc) : start{0}
    {
        extents[0] = nr;
        extents[1] = nc;
        Matrix_impl::compute_strides(*this);
    }

    // Calculate index from a set of subscripts:
    Index operator()(Index i, Index j) const
    {
        return start + i * strides[0] + j * strides[1];
    }

    // Calculate offset given a range.
    template <typename R>
    Index offset(R&& range) const
    {
        constexpr Index zero = 0;
        return start + std::inner_product(strides.begin(), strides.end(),
                                          std::begin(range), zero);
    }

    Index size;                   // total number of elements
    Index start;                  // starting offset
    std::array<Index, 2> extents; // number of elements in each dim
    std::array<Index, 2> strides; // offsets between elements
};

// Matrix_slice to describe three-dimensional matrix (cube).
template <>
struct Matrix_slice<3> {
    // Empty matrix:
    Matrix_slice() = default;

    // Copy semantics:
    Matrix_slice(const Matrix_slice&) = default;
    Matrix_slice& operator=(const Matrix_slice&) = default;

    // Starting offset and extents:
    Matrix_slice(Index s, std::initializer_list<Index> exts) : start(s)
    {
        assert(exts.size() == 3);
        std::copy(exts.begin(), exts.end(), extents.begin());
        Matrix_impl::compute_strides(*this);
    }

    // Starting offset, extents, and strides:
    Matrix_slice(Index s,
                 std::initializer_list<Index> exts,
                 std::initializer_list<Index> strs)
        : start(s)
    {
        assert(exts.size() == 3);
        std::copy(exts.begin(), exts.end(), extents.begin());
        std::copy(strs.begin(), strs.end(), strides.begin());
        size = Matrix_impl::compute_size(extents);
    }

    // N extents:
    Matrix_slice(Index n1, Index n2, Index n3) : start{0}
    {
        extents[0] = n1;
        extents[1] = n2;
        extents[2] = n3;
        Matrix_impl::compute_strides(*this);
    }

    // Calculate index from a set of subscripts:
    Index operator()(Index i, Index j, Index k) const
    {
        return start + i * strides[0] + j * strides[1] + k * strides[2];
    }

    // Calculate offset given a range.
    template <typename R>
    Index offset(R&& range) const
    {
        constexpr Index zero = 0;
        return start + std::inner_product(strides.begin(), strides.end(),
                                          std::begin(range), zero);
    }

    Index size;                   // total number of elements
    Index start;                  // starting offset
    std::array<Index, 3> extents; // number of elements in each dim
    std::array<Index, 3> strides; // offsets between elements
};

//------------------------------------------------------------------------------

// Non-member functions:

// Two Matrix_slices compare equal when they describe the same sequence of
// offsets.
template <std::size_t N>
inline bool operator==(const Matrix_slice<N>& a, const Matrix_slice<N>& b)
{
    return a.start == b.start && a.extents == b.extents &&
           a.strides == b.strides;
}

template <std::size_t N>
inline bool operator!=(const Matrix_slice<N>& a, const Matrix_slice<N>& b)
{
    return !(a == b);
}

//------------------------------------------------------------------------------

// Same extents:

// Return true when two slices describe matrices with the same order and
// extents. The starting offset and strides do not factor into the
// comparison.
//
// An overload is provided for Matrix_type. It compares the descriptors
// of its matrix arguments.

template <std::size_t N>
bool same_extents(const Matrix_slice<N>& a, const Matrix_slice<N>& b)
{
    return a.extents == b.extents;
}

template <typename M1, typename M2>
inline Enable_if<Matrix_type<M1>() && Matrix_type<M2>(), bool>
same_extents(const M1& a, const M2& b)
{
    return same_extents(a.descriptor(), b.descriptor());
}

} // namespace Numlib

#endif // NUMLIB_DENSE_MATRIX_SLICE_H
