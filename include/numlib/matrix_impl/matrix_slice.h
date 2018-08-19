////////////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2018 Stig Rune Sellevag. All rights reserved.
//
// This code is licensed under the MIT License (MIT).
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.
//
////////////////////////////////////////////////////////////////////////////////

#ifndef NUMLIB_MATRIX_MATRIX_SLICE_H
#define NUMLIB_MATRIX_MATRIX_SLICE_H

#include <numlib/matrix_impl/support.h>
#include <numlib/traits/traits.h>
#include <algorithm>
#include <array>
#include <cassert>
#include <initializer_list>


template <std::size_t N>
struct Matrix_slice {
    // Empty matrix:
    Matrix_slice() = default;

    // Copy semantics:
    Matrix_slice(const Matrix_slice&) = default;
    Matrix_slice& operator=(const Matrix_slice&) = default;

    // Starting offset and extents:
    Matrix_slice(std::size_t offset, std::initializer_list<std::size_t> exts);

    // Starting offset, extents, and strides:
    Matrix_slice(std::size_t offset,
                 std::initializer_list<std::size_t> exts,
                 std::initializer_list<std::size_t> strs);

    // N extents:
    template <typename... Dims>
    Matrix_slice(Dims... dims);

    // Calculate index from a set of subscripts:
    // clang-format off
    template <typename... Dims>
#ifdef _MSC_VER // Workaround for internal compiler error in VS 2017
    std::size_t operator()(Dims... dims) const;
#else
    Enable_if<All(Convertible<Dims, std::size_t>()...), std::size_t> 
    operator()(Dims... dims) const;
#endif // _MSC_VER
    // clang-format on

    std::size_t size;                    // total number of elements
    std::size_t start;                   // starting offset
    std::array<std::size_t, N> extents;  // number of elements in each dimension
    std::array<std::size_t, N> strides;  // offsets between elements in each dim
};

template <std::size_t N>
Matrix_slice<N>::Matrix_slice(std::size_t offset,
                              std::initializer_list<std::size_t> exts)
    : start(offset)
{
    assert(exts.size() == N);
    std::copy(exts.begin(), exts.end(), extents.begin());
    Matrix_impl::compute_strides(*this);
}

template <std::size_t N>
Matrix_slice<N>::Matrix_slice(std::size_t offset,
                              std::initializer_list<std::size_t> exts,
                              std::initializer_list<std::size_t> strs)
    : start(offset)
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
    std::size_t args[N]{std::size_t(dims)...};
    std::copy(std::begin(args), std::end(args), extents.begin());
    Matrix_impl::compute_strides(*this);
}

template <std::size_t N>
template <typename... Dims>
#ifdef _MSC_VER  // Workaround for internal compiler error in VS 2017
inline std::size_t
#else
inline Enable_if<All(Convertible<Dims, std::size_t>()...), std::size_t>
#endif  // _MSC_VER
Matrix_slice<N>::operator()(Dims... dims) const
{
    static_assert(sizeof...(Dims) == N,
                  "Matrix_slice<N>::operator(): dimension mismatch");
    std::size_t args[N]{std::size_t(dims)...};
    return start
           + std::inner_product(
                 args, args + N, strides.begin(), std::size_t{0});
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
    Matrix_slice(std::size_t offset, std::size_t exts)
    {
        start      = offset;
        extents[0] = exts;
        strides[0] = 1;
        size       = exts;
    }

    // Starting offset, extents, and strides:
    Matrix_slice(std::size_t offset, std::size_t exts, std::size_t strs)
    {
        start      = offset;
        extents[0] = exts;
        strides[0] = strs;
        size       = exts * strs;
    }

    // N extents:
    Matrix_slice(std::size_t exts)
    {
        start      = 0;
        extents[0] = exts;
        strides[0] = 1;
        size       = exts;
    }

    // Calculate index from a set of subscripts:
    std::size_t operator()(std::size_t i) const { return i; }

    std::size_t size;                    // total number of elements
    std::size_t start;                   // starting offset
    std::array<std::size_t, 1> extents;  // number of elements in each dimension
    std::array<std::size_t, 1> strides;  // offsets between elements in each dim
};

//------------------------------------------------------------------------------

// Non-member functions:

// Return true if the two Matrix_slices have same extents.
template <std::size_t N>
bool same_extents(const Matrix_slice<N>& a, const Matrix_slice<N>& b)
{
    return a.extents == b.extents;
}

#endif  // NUMLIB_MATRIX_MATRIX_SLICE_H
