// Copyright (c) 2008-2010 Kent State University
// Copyright (c) 2011-2012 Texas A&M University
// Copyright (c) 2018 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#ifndef NUMLIB_DENSE_MATRIX_H
#define NUMLIB_DENSE_MATRIX_H

#include <numlib/matrix_impl/matrix_base.h>
#include <initializer_list>
#include <utility>
#include <vector>

namespace Numlib {

// N-dimensional dense matrix class using row-major storage order.
//
// The matrix class provides support for indexing, slicing and basic
// arithmetic operations.
//
// Template parameters:
//   T - The element type stored by the matrix
//   N - The matrix rank
//
template <typename T, std::size_t N>
class Matrix : public Matrix_base<T, N> {
public:
    using size_type = typename Matrix_base<T, N>::size_type;
    using value_type = typename Matrix_base<T, N>::value_type;
    using iterator = typename std::vector<T>::iterator;
    using const_iterator = typename std::vector<T>::const_iterator;

    Matrix() = default;

    // Move construction and assignment:

    Matrix(Matrix&&) = default;
    Matrix& operator=(Matrix&&) = default;

    // Copy construction and assignment:

    Matrix(const Matrix&) = default;
    Matrix& operator=(const Matrix&) = default;

    // Specify the extents.
    template <typename... Exts>
    explicit Matrix(Exts... exts);

    // Construct and assign from Matrix_ref:

    template <typename U>
    Matrix(const Matrix_ref<U, N>&);

    template <typename U>
    Matrix& operator=(const Matrix_ref<U, N>&);

    // Initialize and assign from list:

    template <typename U>
    Matrix(std::initializer_list<U>) = delete; // don't use {} except for elems

    template <typename U>
    Matrix& operator=(std::initializer_list<U>) = delete;

    Matrix(Matrix_initializer<T, N>);
    Matrix& operator=(Matrix_initializer<T, N>);

    ~Matrix() = default;

    // "Flat" element access:
    T* data() { return elems.data(); }
    const T* data() const { return elems.data(); }

    // Properties:

    bool empty() const { return elems.empty(); }

    // Subscripting:

    template <typename... Args>
    Enable_if<Matrix_impl::Requesting_element<Args...>(), T&>
    operator()(Args... args)
    {
        assert(Matrix_impl::check_bounds(this->desc, args...));
        return *(data() + this->desc(args...));
    }

    template <typename... Args>
    Enable_if<Matrix_impl::Requesting_element<Args...>(), const T&>
    operator()(Args... args) const
    {
        assert(Matrix_impl::check_bounds(this->desc, args...));
        return *(data() + this->desc(args...));
    }

    template <typename... Args>
    Enable_if<Matrix_impl::Requesting_slice<Args...>(), Matrix_ref<T, N>>
    operator()(const Args&... args)
    {
        Matrix_slice<N> d;
        d.start = Matrix_impl::do_slice(this->desc, d, args...);
        d.size = Matrix_impl::compute_size(d.extents);
        return {d, data()};
    }

    template <typename... Args>
    Enable_if<Matrix_impl::Requesting_slice<Args...>(), Matrix_ref<const T, N>>
    operator()(const Args&... args) const
    {
        Matrix_slice<N> d;
        d.start = Matrix_impl::do_slice(this->desc, d, args...);
        d.size = Matrix_impl::compute_size(d.extents);
        return {d, data()};
    }

    // Row subscripting.
    Matrix_ref<T, N - 1> operator[](size_type n);
    Matrix_ref<const T, N - 1> operator[](size_type n) const;

    // Return a reference to the n'th row of the matrix.
    Matrix_ref<T, N - 1> row(size_type n);
    Matrix_ref<const T, N - 1> row(size_type n) const;

    // Return a reference to the n'th column of the matrix.
    Matrix_ref<T, N - 1> column(size_type n);
    Matrix_ref<const T, N - 1> column(size_type n) const;

    // Return a reference to the diagonal of a square two-dimensional matrix.
    Matrix_ref<T, N - 1> diag();
    Matrix_ref<const T, N - 1> diag() const;

    // Iterators:

    iterator begin() { return elems.begin(); }
    iterator end() { return elems.end(); }

    const_iterator begin() const { return elems.begin(); }
    const_iterator end() const { return elems.end(); }

    // Mutators:

    void swap(Matrix& m);
    void swap_rows(size_type m, size_type n);

    // Resize matrix (elements not preserved).
    void resize(const Matrix_slice<N>& ms);

    template <typename... Exts>
    void resize(Exts... exts);

    // Apply f(x) for every element x.
    template <typename F>
    Matrix& apply(F f);

    // Apply f(x, mx) for corresponding elements *this and m.
    template <typename M, typename F>
    Enable_if<Matrix_type<M>(), Matrix&> apply(const M& m, F f);

    // Arithmetic operations:

    Matrix& operator=(const T& value); // assignment with scalar

    Matrix& operator+=(const T& value); // scalar addition
    Matrix& operator-=(const T& value); // scalar subtraction
    Matrix& operator*=(const T& value); // scalar multiplication
    Matrix& operator/=(const T& value); // scalar division
    Matrix& operator%=(const T& value); // scalar modulo

    template <typename M>
    Enable_if<Matrix_type<M>(), Matrix&>
    operator+=(const M& m); // matrix addition

    template <typename M>
    Enable_if<Matrix_type<M>(), Matrix&>
    operator-=(const M& m); // matrix subtraction

private:
    std::vector<T> elems;
};

template <typename T, std::size_t N>
template <typename... Exts>
inline Matrix<T, N>::Matrix(Exts... exts)
    : Matrix_base<T, N>{exts...}, elems(this->desc.size)
{
}

template <typename T, std::size_t N>
template <typename U>
inline Matrix<T, N>::Matrix(const Matrix_ref<U, N>& m)
    : Matrix_base<T, N>(m.descriptor()), elems{m.begin(), m.end()}
{
    static_assert(Convertible<U, T>(),
                  "Matrix constructor: incompatible element types");
}

template <typename T, std::size_t N>
template <typename U>
inline Matrix<T, N>& Matrix<T, N>::operator=(const Matrix_ref<U, N>& m)
{
    static_assert(Convertible<U, T>(),
                  "Matrix assignment: incompatible element types");
    this->desc = m.descriptor();
    elems.assign(m.begin(), m.end());
    return *this;
}

template <typename T, std::size_t N>
inline Matrix<T, N>::Matrix(Matrix_initializer<T, N> init)
{
    this->desc.start = 0;
    this->desc.extents = Matrix_impl::derive_extents<N>(init);
    Matrix_impl::compute_strides(this->desc);
    elems.reserve(this->desc.size);
    Matrix_impl::insert_flat(init, elems);
    assert(static_cast<size_type>(elems.size()) == this->desc.size);
}

template <typename T, std::size_t N>
inline Matrix<T, N>& Matrix<T, N>::operator=(Matrix_initializer<T, N> init)
{
    Matrix tmp(init);
    swap(tmp);
    return *this;
}

template <typename T, std::size_t N>
inline Matrix_ref<T, N - 1> Matrix<T, N>::row(size_type n)
{
    assert(n < this->rows());
    auto r = Matrix_impl::slice_dim<0>(this->desc, n);
    return {r, data()};
}

template <typename T, std::size_t N>
inline Matrix_ref<const T, N - 1> Matrix<T, N>::row(size_type n) const
{
    assert(n < this->rows());
    auto r = Matrix_impl::slice_dim<0>(this->desc, n);
    return {r, data()};
}

template <typename T, std::size_t N>
inline Matrix_ref<T, N - 1> Matrix<T, N>::column(size_type n)
{
    assert(n < this->cols());
    auto c = Matrix_impl::slice_dim<1>(this->desc, n);
    return {c, data()};
}

template <typename T, std::size_t N>
inline Matrix_ref<const T, N - 1> Matrix<T, N>::column(size_type n) const
{
    assert(n < this->cols());
    auto c = Matrix_impl::slice_dim<1>(this->desc, n);
    return {c, data()};
}

template <typename T, std::size_t N>
inline Matrix_ref<T, N - 1> Matrix<T, N>::diag()
{
    static_assert(N == 2, "diag: only defined for Matrix of rank 2");
    assert(this->rows() == this->cols());

    Matrix_slice<N - 1> d;
    d.start = this->desc.start;
    d.extents[0] = this->rows();
    d.strides[0] = this->rows() + 1;
    d.size = Matrix_impl::compute_size(d.extents);

    return {d, data()};
}

template <typename T, std::size_t N>
inline Matrix_ref<const T, N - 1> Matrix<T, N>::diag() const
{
    static_assert(N == 2, "diag: only defined for Matrix of rank 2");
    assert(this->rows() == this->cols());

    Matrix_slice<N - 1> d;
    d.start = this->desc.start;
    d.extents[0] = this->rows();
    d.strides[0] = this->rows() + 1;
    d.size = Matrix_impl::compute_size(d.extents);

    return {d, data()};
}

template <typename T, std::size_t N>
inline void Matrix<T, N>::swap(Matrix& m)
{
    std::swap(this->desc, m.desc);
    elems.swap(m.elems);
}

template <typename T, std::size_t N>
inline void Matrix<T, N>::swap_rows(size_type m, size_type n)
{
    auto a = (*this)[m];
    auto b = (*this)[n];
    std::swap_ranges(a.begin(), a.end(), b.begin());
}

template <typename T, std::size_t N>
inline void Matrix<T, N>::resize(const Matrix_slice<N>& ms)
{
    this->desc = ms;
    elems.resize(ms.size);
}

template <typename T, std::size_t N>
template <typename... Exts>
inline void Matrix<T, N>::resize(Exts... exts)
{
    assert(sizeof...(Exts) == this->rank());
    Matrix_slice<N> d{static_cast<size_type>(exts)...}; // avoid C2398 error
    this->desc = d;
    elems.resize(this->size());
}

template <typename T, std::size_t N>
template <typename F>
Matrix<T, N>& Matrix<T, N>::apply(F f)
{
    for (auto& x : elems) {
        f(x);
    }
    return *this;
}

template <typename T, std::size_t N>
template <typename M, typename F>
Enable_if<Matrix_type<M>(), Matrix<T, N>&> Matrix<T, N>::apply(const M& m, F f)
{
    assert(same_extents(this->desc, m.descriptor()));
    auto i = begin();
    auto j = m.begin();
    while (i != end()) {
        f(*i, *j);
        ++i;
        ++j;
    }
    return *this;
}

template <typename T, std::size_t N>
inline Matrix<T, N>& Matrix<T, N>::operator=(const T& value)
{
    return apply([&](T& a) { a = value; });
}

template <typename T, std::size_t N>
inline Matrix<T, N>& Matrix<T, N>::operator+=(const T& value)
{
    return apply([&](T& a) { a += value; });
}

template <typename T, std::size_t N>
inline Matrix<T, N>& Matrix<T, N>::operator-=(const T& value)
{
    return apply([&](T& a) { a -= value; });
}

template <typename T, std::size_t N>
inline Matrix<T, N>& Matrix<T, N>::operator*=(const T& value)
{
    return apply([&](T& a) { a *= value; });
}

template <typename T, std::size_t N>
inline Matrix<T, N>& Matrix<T, N>::operator/=(const T& value)
{
    return apply([&](T& a) { a /= value; });
}

template <typename T, std::size_t N>
inline Matrix<T, N>& Matrix<T, N>::operator%=(const T& value)
{
    return apply([&](T& a) { a %= value; });
}

template <typename T, std::size_t N>
template <typename M>
inline Enable_if<Matrix_type<M>(), Matrix<T, N>&>
// clang-format off
Matrix<T, N>::operator+=(const M& m)
// clang-format on
{
    static_assert(M::order == N, "+=: mismatched Matrix dimensions");
    assert(same_extents(this->desc, m.descriptor()));

    return apply(m, [](T& a, const Value_type<M>& b) { a += b; });
}

template <typename T, std::size_t N>
template <typename M>
inline Enable_if<Matrix_type<M>(), Matrix<T, N>&>
// clang-format off
Matrix<T, N>::operator-=(const M& m)
// clang-format on
{
    static_assert(M::order == N, "-=: mismatched Matrix dimensions");
    assert(same_extents(this->desc, m.descriptor()));

    return apply(m, [](T& a, const Value_type<M>& b) { a -= b; });
}

//------------------------------------------------------------------------------

// Zero-dimensional matrix:
//
// The type Matrix<T, 0> is not really a matrix. It stores a single scalar
// of type T and can only be converted to a reference to that type.
template <typename T>
class Matrix<T, 0> : public Matrix_base<T, 0> {
public:
    Matrix() = default;

    Matrix(const T& x) : Matrix_base<T, 0>(), elem(x) {}

    Matrix& operator=(const T& value)
    {
        elem = value;
        return *this;
    }

    T& operator()() { return elem; }
    const T& operator()() const { return elem; }

    operator T&() { return elem; }
    operator const T&() { return elem; }

private:
    T elem;
};

} // namespace Numlib

#endif // NUMLIB_DENSE_MATRIX_H
