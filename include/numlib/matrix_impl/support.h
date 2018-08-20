// Copyright (c) 2008-2010 Kent State University
// Copyright (c) 2011-2012 Texas A&M University
// Copyright (c) 2018 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#ifndef NUMLIB_MATRIX_SUPPORT_H
#define NUMLIB_MATRIX_SUPPORT_H

#include <algorithm>
#include <array>
#include <cassert>
#include <functional>
#include <numeric>

namespace Numlib {

// Forward declarations:

template <std::size_t N>
struct Matrix_slice;

//------------------------------------------------------------------------------

namespace Matrix_impl {

    template <typename... Args>
    constexpr bool Requesting_element()
    {
        return All(Convertible<Args, std::size_t>()...);
    }

    template <typename... Args>
    constexpr bool Requesting_slice()
    {
        // clang-format off
        return All((Convertible<Args, std::size_t>() || Same<Args, Slice>())...) 
					&& Some(Same<Args, Slice>()...);
        // clang-format on
    }

    //--------------------------------------------------------------------------

    // Matrix list initialization:

    // Forward declaration:

    template <std::size_t N, typename List>
    inline bool check_non_jagged(const List& list);

    template <std::size_t N, typename I, typename List>
    inline Enable_if<(N == 1), void> add_extents(I& first, const List& list)
    {
        *first++ = list.size();
    }

    // Recursion through nested std::initializer_list.
    template <std::size_t N, typename I, typename List>
    inline Enable_if<(N > 1), void> add_extents(I& first, const List& list)
    {
        assert(check_non_jagged<N>(list));
        *first++ = list.size(); // store this size (extent)
        add_extents<N - 1>(first, *list.begin());
    }

    // Determine the shape of the Matrix:
    //   + Checks that the tree is really N deep
    //   + Checks that each row has the same number of elements
    //   + Sets the extent of each row
    //
    template <std::size_t N, typename List>
    inline std::array<std::size_t, N> derive_extents(const List& list)
    {
        std::array<std::size_t, N> a;
        auto f = a.begin();
        add_extents<N>(f, list); // add sizes (extents) to a
        return a;
    }

    // Check that all rows have the same number of elements.
    template <std::size_t N, typename List>
    inline bool check_non_jagged(const List& list)
    {
        auto i = list.begin();
        for (auto j = i + 1; j != list.end(); ++j) {
            if (derive_extents<N - 1>(*i) != derive_extents<N - 1>(*j)) {
                return false;
            }
        }
        return true;
    }

    // When we reach a list with non-initializer_list elements, we insert
    // those elements into our vector.
    template <typename T, typename Vec>
    inline void add_list(const T* first, const T* last, Vec& vec)
    {
        vec.insert(vec.end(), first, last);
    }

    template <typename T, typename Vec>
    inline void add_list(const std::initializer_list<T>* first,
                         const std::initializer_list<T>* last, Vec& vec)
    {
        while (first != last) {
            add_list(first->begin(), first->end(), vec);
            ++first;
        }
    }

    // Copy elements of the tree of std::initializer_list to a Matrix<T, N>.
    template <typename T, typename Vec>
    inline void insert_flat(std::initializer_list<T> list, Vec& vec)
    {
        add_list(list.begin(), list.end(), vec);
    }

    //--------------------------------------------------------------------------

    // Compute strides needed for subscript calculation and the number of
    // elements given the extents.
    //
    // Note: Row-major storage order.
    //
    template <std::size_t N>
    inline void compute_strides(Matrix_slice<N>& ms)
    {
        std::size_t st = 1; // last stride is 1
        for (int i = N - 1; i >= 0; --i) {
            ms.strides[i] = st;
            st *= ms.extents[i];
        }
        ms.size = st;
    }

    // Compute total number of elements given the extents.
    template <std::size_t N>
    inline std::size_t compute_size(const std::array<std::size_t, N>& exts)
    {
        return std::accumulate(exts.begin(), exts.end(), 1,
                               std::multiplies<std::size_t>{});
    }

    // Return true if each element in range is within the bounds of the
    // corresponding extent.
    template <std::size_t N, typename... Dims>
    inline bool check_bounds(const Matrix_slice<N>& slice, Dims... dims)
    {
        std::size_t indexes[N]{std::size_t(dims)...};
        return std::equal(indexes, indexes + N, slice.extents.begin(),
                          std::less<std::size_t>{});
    }

    // Return Matrix_slice describing n'th row.
    template <std::size_t D, std::size_t N>
    Matrix_slice<N - 1> slice_dim(const Matrix_slice<N>& ms, std::size_t n)
    {
        static_assert(N >= 1 && D <= N, "get_row: bad dimension");

        Matrix_slice<N - 1> r;
        r.size = ms.size / ms.extents[D];
        r.start = ms.start + n * ms.strides[D];

        // Copy extents and strides.
        auto i = std::copy_n(ms.extents.begin(), D, r.extents.begin());
        auto j = std::copy_n(ms.strides.begin(), D, r.strides.begin());
        std::copy_n(ms.extents.begin() + D + 1, N - D - 1, i);
        std::copy_n(ms.strides.begin() + D + 1, N - D - 1, j);

        return r;
    }

    // Return starting offset given a slice:

    template <std::size_t D, std::size_t N>
    std::size_t do_slice_dim(const Matrix_slice<N>& os, Matrix_slice<N>& ns,
                             std::size_t s)
    {
        std::size_t i = N - D;
        ns.strides[i] = os.strides[i];
        ns.extents[i] = 1;
        return s * ns.strides[i];
    }

    template <std::size_t D, std::size_t N>
    std::size_t do_slice_dim(const Matrix_slice<N>& os, Matrix_slice<N>& ns,
                             Slice s)
    {
        std::size_t i = N - D;
        ns.strides[i] = s.stride * os.strides[i];
        ns.extents[i] =
            (s.length == std::size_t(-1))
                ? (os.extents[i] - s.start + s.stride - 1) / s.stride
                : s.length;
        return s.start * os.strides[i];
    }

    template <std::size_t N>
    std::size_t do_slice(const Matrix_slice<N>& /* os */,
                         Matrix_slice<N>& /* ns */)
    {
        return 0;
    }

    template <std::size_t N, typename T, typename... Args>
    std::size_t do_slice(const Matrix_slice<N>& os, Matrix_slice<N>& ns,
                         const T& s, const Args&... args)
    {
        std::size_t m = do_slice_dim<sizeof...(Args) + 1>(os, ns, s);
        std::size_t n = do_slice(os, ns, args...);
        return m + n;
    }
} // namespace Matrix_impl

} // namespace Numlib

#endif // NUMLIB_MATRIX_SUPPORT_H
