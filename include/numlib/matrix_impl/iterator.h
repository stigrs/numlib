// Copyright (c) 2008-2010 Kent State University
// Copyright (c) 2011-2012 Texas A&M University
// Copyright (c) 2018 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#ifndef NUMLIB_MATRIX_ITERATOR_H
#define NUMLIB_MATRIX_ITERATOR_H

#include <iterator>
#include <algorithm>

namespace Numlib {

// A slice iterator ranges over the elements of a Matrix_ref specified by
// a Matrix_slice.
//
// Note: Slice_iterator is a forward iterator.
// TODO: Make Slice_iterator bidirectional or random access.
//
template <typename T, std::size_t N>
class Slice_iterator {
public:
    using value_type = Remove_const<T>;
    using reference = T&;
    using pointer = T*;
    using difference_type = std::ptrdiff_t;
    using iterator_category = std::forward_iterator_tag;

    Slice_iterator(const Matrix_slice<N>& ms, T* base, bool limit = false);

private:
    const Matrix_slice<N>& desc; // describes the iterator range
    std::size_t indexes[N];      // counting indexes
    T* ptr;                      // current element
};

template <typename T, std::size_t N>
Slice_iterator<T, N>::Slice_iterator(const Matrix_slice<N>& ms, T* base,
                                     bool limit)
    : desc(ms)
{
    std::fill_n(indexes, N, 0);
    if (limit) {
        indexes[0] = desc.extents[0];
        ptr = base + desc.offset(indexes);
    }
    else {
		ptr = base + ms.start;
	}
}
} // namespace Numlib

#endif // NUMLIB_MATRIX_ITERATOR_H
