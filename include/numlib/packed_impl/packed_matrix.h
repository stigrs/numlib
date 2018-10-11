// Copyright (c) 2018 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#ifndef NUMLIB_PACKED_MATRIX_H
#define NUMLIB_PACKED_MATRIX_H

#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable : 4127) // conditional expression is constant
#endif

#include <numlib/matrix.h>
#include <array>
#include <vector>
#include <algorithm>

namespace Numlib {

// Range-checked packed matrix using row-major storage order.
//
// The packed matrix class provides support for indexing and basic
// arithmetic operations. The storage scheme can be either upper or
// lower triangular.
//
template <typename T, Uplo_scheme Uplo>
class Packed_matrix {
public:
    using value_type = T;
    using size_type = Index;
    using iterator = typename std::vector<T>::iterator;
    using const_iterator = typename std::vector<T>::const_iterator;

    static constexpr Uplo_scheme uplo = Uplo;

    // Empty packed matrix.
    Packed_matrix() = default;

    // Copy semantics:
    Packed_matrix(const Packed_matrix&) = default;
    Packed_matrix& operator=(const Packed_matrix&) = default;

    // Move semantics:
    Packed_matrix(Packed_matrix&&) = default;
    Packed_matrix& operator=(Packed_matrix&&) = default;

    // Construct from extent.
    explicit Packed_matrix(size_type n) : elems(n * (n + 1) / 2), extents{n, n}
    {
    }

    // Construct from C array.
    template <Index np>
    Packed_matrix(size_type n, const T (&ap)[np]);

    // Construct from matrix.
    Packed_matrix(const Matrix<T, 1>& a);
    Packed_matrix(const Matrix<T, 2>& a);

    ~Packed_matrix() = default;

    // "Flat" element access:
    T* data() { return elems.data(); }
    const T* data() const { return elems.data(); }

    // Properties:

    bool empty() const { return elems.empty(); }

    size_type size() const { return elems.size(); }
    size_type rows() const { return extents[0]; }
    size_type cols() const { return extents[1]; }
    size_type extent(size_type dim) const
    {
        assert(0 <= dim && dim < 2);
        return extents[dim];
    }

    // Return UPLO scheme.
    char uplo_scheme() const
    {
        if /* constexpr */ (uplo == upper_triang) { // C++17
            return 'U';
        }
        else { // lower triangular
            return 'L';
        }
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

    void swap(Packed_matrix& ap);
    void resize(size_type n);

    // Apply f(x) for every element x.
    template <typename F>
    Packed_matrix& apply(F f);

    // Arithmetic operations:

    Packed_matrix& operator=(const T& value); // assignment with scalar

    Packed_matrix& operator+=(const T& value); // scalar addition
    Packed_matrix& operator-=(const T& value); // scalar subtraction
    Packed_matrix& operator*=(const T& value); // scalar multiplication
    Packed_matrix& operator/=(const T& value); // scalar division
    Packed_matrix& operator%=(const T& value); // scalar modulo

private:
    std::vector<T> elems;
    std::array<size_type, 2> extents;

    static const T zero;

    T& ref(size_type i, size_type j);
    const T& ref(size_type i, size_type j) const;

    size_type offset(size_type i, size_type j) const;
};

template <typename T, Uplo_scheme Uplo>
template <Index np>
Packed_matrix<T, Uplo>::Packed_matrix(size_type n, const T (&ap)[np])
    : elems(np), extents{n, n}
{
    assert(np >= n * (n + 1) / 2);
    for (size_type i = 0; i < np; ++i) {
        elems[i] = ap[i];
    }
}

template <typename T, Uplo_scheme Uplo>
Packed_matrix<T, Uplo>::Packed_matrix(const Matrix<T, 1>& a) : elems(a.size())
{
    size_type n = narrow_cast<size_type>(
        std::ceil(-1.0 + std::sqrt(1.0 + 4.0 * 2.0 * a.size()) / 2.0));
    assert(a.size() >= n * (n + 1) / 2);
    extents[0] = n;
    extents[1] = n;
    for (size_type i = 0; i < size(); ++i) {
        elems[i] = a(i);
    }
}

template <typename T, Uplo_scheme Uplo>
Packed_matrix<T, Uplo>::Packed_matrix(const Matrix<T, 2>& a)
    : elems(a.rows() * (a.rows() + 1) / 2), extents{a.rows(), a.cols()}
{
    assert(a.rows() == a.cols());
    if /* constexpr */ (Uplo == upper_triang) { // C++17
        for (size_type i = 0; i < a.rows(); ++i) {
            for (size_type j = i; j < a.cols(); ++j) {
                (*this)(i, j) = a(i, j);
            }
        }
    }
    else { // lower triangular
        for (size_type i = 0; i < a.rows(); ++i) {
            for (size_type j = 0; j <= i; ++j) {
                (*this)(i, j) = a(i, j);
            }
        }
    }
}

template <typename T, Uplo_scheme Uplo>
inline void Packed_matrix<T, Uplo>::swap(Packed_matrix<T, Uplo>& ap)
{
    static_assert(uplo == ap.uplo,
                  "packed matrices have different triangular form");

    elems.swap(ap.elems);
    std::swap(extents, ap.extents);
}

template <typename T, Uplo_scheme Uplo>
inline void Packed_matrix<T, Uplo>::resize(size_type n)
{
    elems.resize(n * (n + 1) / 2);
    extents = {n, n};
}

template <typename T, Uplo_scheme Uplo>
template <typename F>
Packed_matrix<T, Uplo>& Packed_matrix<T, Uplo>::apply(F f)
{
    for (auto& x : elems) {
        f(x);
    }
    return *this;
}

template <typename T, Uplo_scheme Uplo>
inline Packed_matrix<T, Uplo>& Packed_matrix<T, Uplo>::operator=(const T& value)
{
    return apply([&](T& a) { a = value; });
}

template <typename T, Uplo_scheme Uplo>
inline Packed_matrix<T, Uplo>& Packed_matrix<T, Uplo>::
operator+=(const T& value)
{
    return apply([&](T& a) { a += value; });
}

template <typename T, Uplo_scheme Uplo>
inline Packed_matrix<T, Uplo>& Packed_matrix<T, Uplo>::
operator-=(const T& value)
{
    return apply([&](T& a) { a -= value; });
}

template <typename T, Uplo_scheme Uplo>
inline Packed_matrix<T, Uplo>& Packed_matrix<T, Uplo>::
operator*=(const T& value)
{
    return apply([&](T& a) { a *= value; });
}

template <typename T, Uplo_scheme Uplo>
inline Packed_matrix<T, Uplo>& Packed_matrix<T, Uplo>::
operator/=(const T& value)
{
    return apply([&](T& a) { a /= value; });
}

template <typename T, Uplo_scheme Uplo>
inline Packed_matrix<T, Uplo>& Packed_matrix<T, Uplo>::
operator%=(const T& value)
{
    return apply([&](T& a) { a %= value; });
}

template <typename T, Uplo_scheme Uplo>
inline T& Packed_matrix<T, Uplo>::ref(size_type i, size_type j)
{
    static_assert(Uplo == upper_triang || Uplo == lower_triang,
                  "Packed_matrix: bad storage scheme");

    assert(0 <= i && i < extents[0]);
    assert(0 <= j && j < extents[1]);
    if /* constexpr */ (Uplo == upper_triang) { // C++17
        assert(i <= j);
        return elems[offset(i, j)];
    }
    else if /* constexpr */ (Uplo == lower_triang) { // C++17
        assert(j <= i);
        return elems[offset(i, j)];
    }
}

template <typename T, Uplo_scheme Uplo>
inline const T& Packed_matrix<T, Uplo>::ref(size_type i, size_type j) const
{
    static_assert(Uplo == upper_triang || Uplo == lower_triang,
                  "Packed_matrix: bad storage scheme");

    assert(0 <= i && i < extents[0]);
    assert(0 <= j && j < extents[1]);
    if /* constexpr */ (Uplo == upper_triang) { // C++17
        if (i <= j) {
            return elems[offset(i, j)];
        }
    }
    else if /* constexpr */ (Uplo == lower_triang) { // C++17
        if (j <= i) {
            return elems[offset(i, j)];
        }
    }
    return zero;
}

template <typename T, Uplo_scheme Uplo>
inline Index Packed_matrix<T, Uplo>::offset(size_type i, size_type j) const
{
    static_assert(Uplo == lower_triang || Uplo == upper_triang,
                  "Packed_matrix: bad storage scheme");

    size_type res = 0;
    if /* constexpr */ (Uplo == upper_triang) { // C++17
        res = j + i * (2 * extents[0] - i - 1) / 2;
    }
    else if /* constexpr */ (Uplo == lower_triang) { // C++17
        res = j + i * (i + 1) / 2;
    }
    return res;
}

// clang-format off
template <typename T, Uplo_scheme Uplo>
const typename Packed_matrix<T, Uplo>::value_type 
Packed_matrix<T, Uplo>::zero = value_type{0};
// clang-format on

} // namespace Numlib

#ifdef _MSC_VER
#pragma warning(pop)
#endif

#endif // NUMLIB_PACKED_MATRIX_H
