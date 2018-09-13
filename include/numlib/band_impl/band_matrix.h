// Copyright (c) 2018 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#ifndef NUMLIB_BAND_MATRIX_BAND_MATRIX_H
#define NUMLIB_BAND_MATRIX_BAND_MATRIX_H

#include <array>
#include <vector>
#include <algorithm>
#include <numlib/matrix.h>

namespace Numlib {

// Range-checked band matrix using column-major storage order.
//
// The band matrix class provides support for indexing and basic
// arithmetic operations.
//
// Note:
// Column-major storage order was selected in order to enable use of
// Intel MKL routines.
//
// Template paramenters:
//   T - The element type stored by the band matrix
//
template <typename T>
class Band_matrix {
public:
    using value_type = T;
    using size_type = Index;
    using iterator = typename std::vector<T>::iterator;
    using const_iterator = typename std::vector<T>::const_iterator;

    // Empty band matrix.
    Band_matrix() = default;

    // Copy semantics:
    Band_matrix(const Band_matrix&) = default;
    Band_matrix& operator=(const Band_matrix&) = default;

    // Move semantics:
    Band_matrix(Band_matrix&&) = default;
    Band_matrix& operator=(Band_matrix&&) = default;

    // Construct from extents.
    Band_matrix(size_type m, size_type n, size_type kl, size_type ku)
        : elems((kl + ku + 1) * n), extents{m, n}, bwidth{kl, ku}
    {
    }

    // Construct from C array.
    template <Index nb>
    Band_matrix(size_type m,
                size_type n,
                size_type kl,
                size_type ku,
                const T (&ab)[nb]);

    // Construct from dense matrix.
    Band_matrix(size_type kl, size_type ku, const Matrix<T, 2>& a);

    // Destructor.
    ~Band_matrix() = default;

    // "Flat" element access:
    T* data() { return elems.data(); }
    const T* data() const { return elems.data(); }

    // Properties:

    bool empty() const { return elems.empty(); }

    size_type size() const { return elems.size(); }
    size_type rows() const { return extents[0]; }
    size_type cols() const { return extents[1]; }
    size_type leading_dim() const { return bwidth[0] + bwidth[1] + 1; }
    size_type lower() const { return bwidth[0]; }
    size_type upper() const { return bwidth[1]; }

    size_type extent(size_type dim) const
    {
        assert(dim >= 0 && dim < 2);
        return extents[dim];
    }

    // Subscripting:

    T& operator()(size_type i, size_type j) { return ref(i, j); }
    const T& operator()(size_type i, size_type j) const { return ref(i, j); }

    // Iterators:

    iterator begin() { return elems.begin(); }
    const_iterator begin() const { return elems.begin(); }

    iterator end() { return elems.end(); }
    const_iterator end() const { return elems.end(); }

    // Mutators:

    void swap(Band_matrix& ab);
    void resize(size_type m, size_type n, size_type kl, size_type ku);

    // Apply f(x) for every element x.
    template <typename F>
    Band_matrix& apply(F f);

    // Arithmetic operations:

    Band_matrix& operator=(const T& value); // assignment with scalar

    Band_matrix& operator+=(const T& value); // scalar addition
    Band_matrix& operator-=(const T& value); // scalar subtraction
    Band_matrix& operator*=(const T& value); // scalar multiplication
    Band_matrix& operator/=(const T& value); // scalar division
    Band_matrix& operator%=(const T& value); // scalar modulo

private:
    std::vector<T> elems;
    std::array<size_type, 2> extents;
    std::array<size_type, 2> bwidth;

    static const T zero;

    T& ref(size_type i, size_type j);
    const T& ref(size_type i, size_type j) const;

    size_type offset(size_type i, size_type j) const;
};

template <typename T>
template <Index nb>
Band_matrix<T>::Band_matrix(
    size_type m, size_type n, size_type kl, size_type ku, const T (&ab)[nb])
    : elems(nb), extents{m, n}, bwidth{kl, ku}
{
    assert(nb >= (kl + ku + 1) * n);
    for (size_type i = 0; i < nb; ++i) {
        elems[i] = ab[i];
    }
}

template <typename T>
Band_matrix<T>::Band_matrix(size_type kl, size_type ku, const Matrix<T, 2>& a)
    : elems((kl + ku + 1) * a.cols()),
      extents{a.rows(), a.cols()},
      bwidth{kl, ku}
{
    for (size_type j = 0; j < a.cols(); ++j) {
        for (size_type i = std::max(size_type{0}, j - ku);
             i < std::min(a.rows(), j + kl + 1); ++i) {
            elems[offset(i, j)] = a(i, j);
        }
    }
}

template <typename T>
inline void Band_matrix<T>::swap(Band_matrix<T>& ab)
{
    elems.swap(ab.elems);
    std::swap(extents, ab.extents);
    std::swap(bwidth, ab.bwidth);
}

template <typename T>
inline void
Band_matrix<T>::resize(size_type m, size_type n, size_type kl, size_type ku)
{
    elems.resize((kl + ku + 1) * n);
    extents = {m, n};
    bwidth = {kl, ku};
}

template <typename T>
template <typename F>
Band_matrix<T>& Band_matrix<T>::apply(F f)
{
    for (auto& x : elems) {
        f(x);
    }
    return *this;
}

template <typename T>
inline Band_matrix<T>& Band_matrix<T>::operator=(const T& value)
{
    return apply([&](T& a) { a = value; });
}

template <typename T>
inline Band_matrix<T>& Band_matrix<T>::operator+=(const T& value)
{
    return apply([&](T& a) { a += value; });
}

template <typename T>
inline Band_matrix<T>& Band_matrix<T>::operator-=(const T& value)
{
    return apply([&](T& a) { a -= value; });
}

template <typename T>
inline Band_matrix<T>& Band_matrix<T>::operator*=(const T& value)
{
    return apply([&](T& a) { a *= value; });
}

template <typename T>
inline Band_matrix<T>& Band_matrix<T>::operator/=(const T& value)
{
    return apply([&](T& a) { a /= value; });
}

template <typename T>
inline Band_matrix<T>& Band_matrix<T>::operator%=(const T& value)
{
    return apply([&](T& a) { a %= value; });
}

template <typename T>
inline T& Band_matrix<T>::ref(size_type i, size_type j)
{
    assert(i >= 0 && i < extents[0]);
    assert(j >= 0 && j < extents[1]);

    assert(std::max(size_type{0}, j - bwidth[1]) <= i &&
           i < std::min(extents[1], j + bwidth[0] + 1));
    return elems[offset(i, j)];
}

template <typename T>
inline const T& Band_matrix<T>::ref(size_type i, size_type j) const
{
    assert(i >= 0 && i < extents[0]);
    assert(j >= 0 && j < extents[1]);

    if (std::max(size_type{0}, j - bwidth[1]) <= i &&
        i < std::min(extents[1], j + bwidth[0] + 1)) {
        return elems[offset(i, j)];
    }
    else {
        return zero;
    }
}

template <typename T>
inline Index Band_matrix<T>::offset(size_type i, size_type j) const
{
    return bwidth[1] + i - j + j * leading_dim();
}

template <typename T>
const typename Band_matrix<T>::value_type Band_matrix<T>::zero = value_type{0};

} // namespace Numlib

#endif // NUMLIB_BAND_MATRIX_BAND_MATRIX_H
