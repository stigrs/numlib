// Copyright (c) 2018 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#ifndef NUMLIB_SPARSE_VECTOR_H
#define NUMLIB_SPARSE_VECTOR_H

#include <vector>
#include <initializer_list>

namespace Numlib {

template <typename T>
class Sparse_vector {
public:
    using value_type = T;
    using size_type = std::ptrdiff_t;
    using iterator = typename std::vector<T>::iterator;
    using const_iterator = typename std::vector<T>::const_iterator;

    // Empty sparse vector.
    Sparse_vector() = default;

    // Copy semantics:
    Sparse_vector(const Sparse_vector&) = default;
    Sparse_vector& operator=(const Sparse_vector&) = default;

    // Move semantics:
    Sparse_vector(Sparse_vector&&) = default;
    Sparse_vector& operator=(Sparse_vector&&) = default;

    // Construct from extents.
    explicit Sparse_vector(size_type n) : elems(n), indx(n) {}

    // Construct from values and indexes:

    Sparse_vector(const std::vector<T>& val, const std::vector<size_type>& loc)
        : elems(val), indx(loc)
    {
    }

    template <std::ptrdiff_t n>
    Sparse_vector(const T (&val)[n], const std::ptrdiff_t (&loc)[n]);

    ~Sparse_vector() = default;

private:
    std::vector<T> elems;
    std::vector<size_type> indx;
    static const T zero;
};

template <typename T>
template <std::ptrdiff_t n>
Sparse_vector<T>::Sparse_vector(const T (&val)[n],
                                const std::ptrdiff_t (&loc)[n])
    : elems(n), indx(n)
{
    for (size_type i = 0; i < n; ++i) {
        elems[i] = val[i];
        indx[i] = loc[i];
    }
}

// clang-format off
template <typename T>
const typename Sparse_vector<T>::value_type 
Sparse_vector<T>::zero = value_type{0};
// clang-format on

} // namespace Numlib

#endif // NUMLIB_SPARSE_VECTOR_H
