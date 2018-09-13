// Copyright (c) 2008-2010 Kent State University
// Copyright (c) 2011-2012 Texas A&M University
// Copyright (c) 2018 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#ifndef NUMLIB_DENSE_MATRIX_ITERATOR_H
#define NUMLIB_DENSE_MATRIX_ITERATOR_H

#include <iterator>
#include <algorithm>
#include <array>

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

    Slice_iterator(const Matrix_slice<N>& ms, T* p, bool limit = false);
    Slice_iterator& operator=(const Slice_iterator& it);

    // Return iterators describing the slice.
    const Matrix_slice<N>& descriptor() const { return desc; }

    // Readable:

    T& operator*() { return *ptr; }
    T* operator->() { return ptr; }

    const T& operator*() const { return *ptr; }
    const T* operator->() const { return ptr; }

    // Forward iterator:

    Slice_iterator& operator++()
    {
        increment();
        return *this;
    }

    Slice_iterator operator++(int)
    {
        Slice_iterator s = *this;
        increment();
        return s;
    }

private:
    // Move to the next element in the range.
    void increment();

    const Matrix_slice<N>& desc;  // describes the iterator range
    std::array<Index, N> indexes; // counting indexes
    T* ptr;                       // current element
};

template <typename T, std::size_t N>
Slice_iterator<T, N>::Slice_iterator(const Matrix_slice<N>& ms,
                                     T* p,
                                     bool limit)
    : desc(ms)
{
    std::fill(indexes.begin(), indexes.end(), 0);
    if (limit) {
        indexes[0] = desc.extents[0];
        ptr = p + desc.offset(indexes);
    }
    else {
        ptr = p + ms.start;
    }
}

template <typename T, std::size_t N>
Slice_iterator<T, N>& Slice_iterator<T, N>::operator=(const Slice_iterator& it)
{
    std::copy(it.indexes.begin(), it.indexes.end(), indexes.begin());
    ptr = it.ptr;
    return *this;
}

template <typename T, std::size_t N>
void Slice_iterator<T, N>::increment()
{
    Index d = N - 1;
    while (true) {
        ptr += desc.strides[d];
        ++indexes[d];

        // If have not yet counted to the extent of the current dimension, then
        // continue to do so in the nex iteration.
        if (indexes[d] != desc.extents[d]) {
            break;
        }

        // Otherwise, if we have not counted to the extent in the outermost
        // dimension, move to the next dimension and try again. If d == 0, then
        // we have counted through the entire slice.
        if (d != 0) {
            ptr -= desc.strides[d] * desc.extents[d];
            indexes[d] = 0;
            --d;
        }
        else {
            break;
        }
    }
}

//------------------------------------------------------------------------------

// Equality comparable:

// Two slice iterators are equality comparable when their slices compute the
// same sequence of elements and the iterators refer to the same element.
template <typename T, std::size_t N>
inline bool operator==(const Slice_iterator<T, N>& a,
                       const Slice_iterator<T, N>& b)
{
    assert(a.descriptor() == b.descriptor());
    return &*a == &*b;
}

template <typename T, std::size_t N>
inline bool operator!=(const Slice_iterator<T, N>& a,
                       const Slice_iterator<T, N>& b)
{
    return !(a == b);
}

} // namespace Numlib

#endif // NUMLIB_DENSE_MATRIX_ITERATOR_H
