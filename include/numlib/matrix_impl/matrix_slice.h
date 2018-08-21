// Copyright (c) 2008-2010 Kent State University
// Copyright (c) 2011-2012 Texas A&M University
// Copyright (c) 2018 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#ifndef NUMLIB_MATRIX_MATRIX_SLICE_H
#define NUMLIB_MATRIX_MATRIX_SLICE_H

#include <algorithm>
#include <array>
#include <cassert>
#include <initializer_list>
#include <iterator>

namespace num {

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
    Matrix_slice() = default;

    // Copy semantics:
    Matrix_slice(const Matrix_slice&) = default;
    Matrix_slice& operator=(const Matrix_slice&) = default;

    // Starting offset and extents:
    Matrix_slice(std::size_t s, std::initializer_list<std::size_t> exts);

    // Starting offset, extents, and strides:
    Matrix_slice(std::size_t s, std::initializer_list<std::size_t> exts,
                 std::initializer_list<std::size_t> strs);

    // N extents:
    template <typename... Dims>
    Matrix_slice(Dims... dims);

    // Calculate index from a set of subscripts:
    template <typename... Dims>
#ifdef _MSC_VER // Workaround for internal compiler error in VS 2017
    std::size_t operator()(Dims... dims) const;
#else
    Enable_if<All(Convertible<Dims, std::size_t>()...), std::size_t>
    operator()(Dims... dims) const;
#endif // _MSC_VER

    // Calculate offset given a range.
    template <typename R>
    std::size_t offset(R&& range) const;

    std::size_t size;                   // total number of elements
    std::size_t start;                  // starting offset
    std::array<std::size_t, N> extents; // number of elements in each dimension
    std::array<std::size_t, N> strides; // offsets between elements in each dim
};

template <std::size_t N>
Matrix_slice<N>::Matrix_slice(std::size_t s,
                              std::initializer_list<std::size_t> exts)
    : start(s)
{
    assert(exts.size() == N);
    std::copy(exts.begin(), exts.end(), extents.begin());
    matrix_impl::compute_strides(*this);
}

template <std::size_t N>
Matrix_slice<N>::Matrix_slice(std::size_t s,
                              std::initializer_list<std::size_t> exts,
                              std::initializer_list<std::size_t> strs)
    : start(s)
{
    assert(exts.size() == N);
    std::copy(exts.begin(), exts.end(), extents.begin());
    std::copy(strs.begin(), strs.end(), strides.begin());
    size = matrix_impl::compute_size(extents);
}

template <std::size_t N>
template <typename... Dims>
Matrix_slice<N>::Matrix_slice(Dims... dims) : start{0}
{
    static_assert(sizeof...(Dims) == N,
                  "Matrix_slice<N>::Matrix_slice(Dims...): dimension mismatch");
    std::size_t args[N]{std::size_t(dims)...};
    std::copy(std::begin(args), std::end(args), extents.begin());
    matrix_impl::compute_strides(*this);
}

template <std::size_t N>
template <typename... Dims>
#ifdef _MSC_VER // Workaround for internal compiler error in VS 2017
inline std::size_t
#else
inline Enable_if<All(Convertible<Dims, std::size_t>()...), std::size_t>
#endif // _MSC_VER
Matrix_slice<N>::operator()(Dims... dims) const
{
    static_assert(sizeof...(Dims) == N,
                  "Matrix_slice<N>::operator(): dimension mismatch");
    std::size_t args[N]{std::size_t(dims)...};
    return start +
           std::inner_product(args, args + N, strides.begin(), std::size_t{0});
}

template <std::size_t N>
template <typename R>
inline std::size_t Matrix_slice<N>::offset(R&& range) const
{
    constexpr std::size_t zero = 0;
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
    Matrix_slice(std::size_t s, std::size_t exts)
    {
        start = s;
        extents[0] = exts;
        strides[0] = 1;
        size = exts;
    }

    // Starting offset, extents, and strides:
    Matrix_slice(std::size_t s, std::size_t exts, std::size_t strs)
    {
        start = s;
        extents[0] = exts;
        strides[0] = strs;
        size = exts * strs;
    }

    // N extents:
    Matrix_slice(std::size_t exts)
    {
        start = 0;
        extents[0] = exts;
        strides[0] = 1;
        size = exts;
    }

    // Calculate index from a set of subscripts:
    std::size_t operator()(std::size_t i) const
    {
        return start + i * strides[0];
    }

    // Calculate offset given a range.
    template <typename R>
    std::size_t offset(R&& range) const
    {
        constexpr std::size_t zero = 0;
        return start + std::inner_product(strides.begin(), strides.end(),
                                          std::begin(range), zero);
    }

    std::size_t size;                   // total number of elements
    std::size_t start;                  // starting offset
    std::array<std::size_t, 1> extents; // number of elements in each dimension
    std::array<std::size_t, 1> strides; // offsets between elements in each dim
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
    Matrix_slice(std::size_t s, std::initializer_list<std::size_t> exts)
        : start(s)
    {
        assert(exts.size() == 2);
        std::copy(exts.begin(), exts.end(), extents.begin());
        matrix_impl::compute_strides(*this);
    }

    // Starting offset, extents, and strides:
    Matrix_slice(std::size_t s, std::initializer_list<std::size_t> exts,
                 std::initializer_list<std::size_t> strs)
        : start(s)
    {
        assert(exts.size() == 2);
        std::copy(exts.begin(), exts.end(), extents.begin());
        std::copy(strs.begin(), strs.end(), strides.begin());
        size = matrix_impl::compute_size(extents);
    }

    // N extents:
    Matrix_slice(std::size_t nr, std::size_t nc) : start{0}
    {
        extents[0] = nr;
        extents[1] = nc;
        matrix_impl::compute_strides(*this);
    }

    // Calculate index from a set of subscripts:
    std::size_t operator()(std::size_t i, std::size_t j) const
    {
        return start + i * strides[0] + j * strides[1];
    }

    // Calculate offset given a range.
    template <typename R>
    std::size_t offset(R&& range) const
    {
        constexpr std::size_t zero = 0;
        return start + std::inner_product(strides.begin(), strides.end(),
                                          std::begin(range), zero);
    }

    std::size_t size;                   // total number of elements
    std::size_t start;                  // starting offset
    std::array<std::size_t, 2> extents; // number of elements in each dimension
    std::array<std::size_t, 2> strides; // offsets between elements in each dim
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
    Matrix_slice(std::size_t s, std::initializer_list<std::size_t> exts)
        : start(s)
    {
        assert(exts.size() == 3);
        std::copy(exts.begin(), exts.end(), extents.begin());
        matrix_impl::compute_strides(*this);
    }

    // Starting offset, extents, and strides:
    Matrix_slice(std::size_t s, std::initializer_list<std::size_t> exts,
                 std::initializer_list<std::size_t> strs)
        : start(s)
    {
        assert(exts.size() == 3);
        std::copy(exts.begin(), exts.end(), extents.begin());
        std::copy(strs.begin(), strs.end(), strides.begin());
        size = matrix_impl::compute_size(extents);
    }

    // N extents:
    Matrix_slice(std::size_t n1, std::size_t n2, std::size_t n3) : start{0}
    {
        extents[0] = n1;
        extents[1] = n2;
        extents[2] = n3;
        matrix_impl::compute_strides(*this);
    }

    // Calculate index from a set of subscripts:
    std::size_t operator()(std::size_t i, std::size_t j, std::size_t k) const
    {
        return start + i * strides[0] + j * strides[1] + k * strides[2];
    }

    // Calculate offset given a range.
    template <typename R>
    std::size_t offset(R&& range) const
    {
        constexpr std::size_t zero = 0;
        return start + std::inner_product(strides.begin(), strides.end(),
                                          std::begin(range), zero);
    }

    std::size_t size;                   // total number of elements
    std::size_t start;                  // starting offset
    std::array<std::size_t, 3> extents; // number of elements in each dimension
    std::array<std::size_t, 3> strides; // offsets between elements in each dim
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

// Return true when two slices describe matrices with the same rank and
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

} // namespace num

#endif // NUMLIB_MATRIX_MATRIX_SLICE_H
