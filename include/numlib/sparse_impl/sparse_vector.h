// Copyright (c) 2018 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#ifndef NUMLIB_SPARSE_VECTOR_H
#define NUMLIB_SPARSE_VECTOR_H

#include <algorithm>
#include <vector>
#include <initializer_list>
#include <numlib/matrix.h>

namespace Numlib {

// Range-checked sparse vector class.
//
// This class provides a basic framework for implementing sparse vector
// methods that utilize the Intel Math Kernel Library.
//
// Note:
// - It is assumed that the sparse vector is initialized with elements indices
//   sorted in ascending order.
// - Zero and one-based indexing are supported.
// - New elements are inserted so that the index order is preserved.
//
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

    // Construct from values and indexes:

    Sparse_vector(const std::vector<T>& val, const std::vector<size_type>& loc)
        : elems(val), indx(loc)
    {
    }

    template <std::ptrdiff_t n>
    Sparse_vector(const T (&val)[n], const std::ptrdiff_t (&loc)[n]);

    // Construct from initializer list.
    Sparse_vector(std::initializer_list<std::pair<size_type, T>> list);

    // Assign from values and indexes.
    Sparse_vector&
    operator=(std::initializer_list<std::pair<size_type, T>> list);

    ~Sparse_vector() = default;

    // Flat" element access:

    T* data() { return elems.data(); }
    const T* data() const { return elems.data(); }

    // Access underlying arrays:

    auto& values() { return elems; }
    const auto& values() const { return elems; }
    const auto& index() const { return indx; }
    const auto& index_zero_based() const { return indx; }
    auto index_one_based() const;

    // Properties:

    bool empty() const { return elems.empty(); }

    size_type num_nonzero() const { return elems.size(); }
    size_type size() const
    {
        return *std::max_element(indx.begin(), indx.end()) + 1;
    }

    // Subscripting:

    size_type loc(size_type i) const
    {
        assert(0 <= i && i < num_nonzero());
        return indx[i];
    }

    const T& operator()(size_type i) const { return ref(i); }

    // Iterators:

    iterator begin() { return elems.begin(); }
    const_iterator begin() const { return elems.begin(); }

    iterator end() { return elems.end(); }
    const_iterator end() const { return elems.end(); }

    // Mutators:

    void swap(Sparse_vector& v);
    void insert(size_type i, const T& val);

    // Apply f(x) for every element x.
    template <typename F>
    Sparse_vector& apply(F f);

    // Arithmetic operations:

    Sparse_vector& operator*=(const T& value); // scalar multiplication
    Sparse_vector& operator/=(const T& value); // scalar division

private:
    std::vector<T> elems;
    std::vector<size_type> indx;
    static const T zero;

    const T& ref(size_type i) const;
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

template <typename T>
Sparse_vector<T>::Sparse_vector(
    std::initializer_list<std::pair<size_type, T>> list)
    : elems(list.size()), indx(list.size())
{
    size_type i = 0;
    for (const auto& il : list) {
        elems[i] = std::get<1>(il);
        indx[i] = std::get<0>(il);
        ++i;
    }
}

template <typename T>
Sparse_vector<T>& Sparse_vector<T>::
operator=(std::initializer_list<std::pair<size_type, T>> list)
{
    elems.resize(list.size());
    indx.resize(list.size());
    size_type i = 0;
    for (const auto& il : list) {
        elems[i] = std::get<1>(il);
        indx[i] = std::get<0>(il);
        ++i;
    }
    return *this;
}

template <typename T>
auto Sparse_vector<T>::index_one_based() const
{
    auto res = indx;
    for (auto& i : res) {
        i += 1;
    }
    return res;
}

template <typename T>
inline void Sparse_vector<T>::swap(Sparse_vector<T>& v)
{
    elems.swap(v.elems);
    indx.swap(v.indx);
}

template <typename T>
inline void Sparse_vector<T>::insert(size_type i, const T& val)
{
    if (val != T{0}) { // zero values should not be stored
        if (std::find(indx.begin(), indx.end(), i) == indx.end()) {
            // do not replace any existing elements
            auto pos = std::upper_bound(indx.begin(), indx.end(), i);
            size_type offset = std::distance(indx.begin(), pos);
            elems.insert(elems.begin() + offset, val);
            indx.insert(pos, i);
        }
    }
}

template <typename T>
template <typename F>
inline Sparse_vector<T>& Sparse_vector<T>::apply(F f)
{
    for (auto& x : elems) {
        f(x);
    }
    return *this;
}

template <typename T>
inline Sparse_vector<T>& Sparse_vector<T>::operator*=(const T& value)
{
    return apply([&](T& a) { a *= value; });
}

template <typename T>
inline Sparse_vector<T>& Sparse_vector<T>::operator/=(const T& value)
{
    return apply([&](T& a) { a /= value; });
}

template <typename T>
inline const T& Sparse_vector<T>::ref(size_type i) const
{
    assert(0 <= i && i < size());

    auto pos = std::find(indx.begin(), indx.end(), i);
    size_type offset = std::distance(indx.begin(), pos);
    return offset >= 0 && offset < num_nonzero() ? elems[offset] : zero;
}

// clang-format off
template <typename T>
const typename Sparse_vector<T>::value_type 
Sparse_vector<T>::zero = value_type{0};
// clang-format on

} // namespace Numlib

#endif // NUMLIB_SPARSE_VECTOR_H
