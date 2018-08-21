// Copyright (c) 2008-2010 Kent State University
// Copyright (c) 2011-2012 Texas A&M University
// Copyright (c) 2018 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#ifndef NUMLIB_MATRIX_MATRIX_REF_H
#define NUMLIB_MATRIX_MATRIX_REF_H

#include <numlib/matrix_impl/matrix_base.h>
#include <utility>

namespace num {

//------------------------------------------------------------------------------

template <typename T, std::size_t N>
class Matrix_ref : public Matrix_base<T, N> {
public:
    using size_type = typename Matrix_base<T, N>::size_type;
    using value_type = Remove_const<T>;
    using iterator = Slice_iterator<T, N>;
    using const_iterator = Slice_iterator<const T, N>;

    Matrix_ref() = default;

    // Move construction and assignment:
    Matrix_ref(Matrix_ref&&);
    Matrix_ref& operator=(Matrix_ref&&);

    // Copy construction and assignment:
    Matrix_ref(const Matrix_ref&);
    Matrix_ref& operator=(const Matrix_ref&);

    // Slice initialization:
    //
    // Initialize the Matrix_ref over the Matrix_slice, starting at the
    // element pointed to by p.
    Matrix_ref(const Matrix_slice<N>& ms, T* p) : Matrix_base<T, N>(ms), ptr(p)
    {
    }

    // Matrix initialization:
    //
    // Initialize the Matrix_ref so that it refers to another Matrix m, or
    // assign the elements of that Matrix into this sub-matrix.
    Matrix_ref(Matrix<value_type, N>& m);
    Matrix_ref(const Matrix<value_type, N>& m);
    Matrix_ref(Matrix<value_type, N>&&) = delete; // avoid memory leaks

    Matrix_ref& operator=(const Matrix<value_type, N>& m);

    // Sub-matrix conversion:
    //
    // Allow implict conversion from a non-const Matrix_ref to a const
    // Matrix_ref.
    template <typename U>
    Matrix_ref(const Matrix_ref<U, N>& m);

    template <typename U>
    Matrix_ref& operator=(const Matrix_ref<U, N>& m);

    ~Matrix_ref() = default;

    // "Flat" element access:
    T* data() { return ptr; }
    const T* data() const { return ptr; }

    // Subscripting:

    template <typename... Args>
    Enable_if<matrix_impl::Requesting_element<Args...>(), T&>
    operator()(Args... args)
    {
        assert(matrix_impl::check_bounds(this->desc, args...));
        return *(data() + this->desc(args...));
    }

    template <typename... Args>
    Enable_if<matrix_impl::Requesting_element<Args...>(), const T&>
    operator()(Args... args) const
    {
        assert(matrix_impl::check_bounds(this->desc, args...));
        return *(data() + this->desc(args...));
    }

    template <typename... Args>
    Enable_if<matrix_impl::Requesting_slice<Args...>(), Matrix_ref<T, N>>
    operator()(Args... args)
    {
        Matrix_slice<N> d;
        d.start = matrix_impl::do_slice(this->desc, d, args...);
        d.size = matrix_impl::compute_size(d.extents);
        return {d, data()};
    }

    template <typename... Args>
    Enable_if<matrix_impl::Requesting_slice<Args...>(), Matrix_ref<const T, N>>
    operator()(Args... args) const
    {
        Matrix_slice<N> d;
        d.start = matrix_impl::do_slice(this->desc, d, args...);
        d.size = matrix_impl::compute_size(d.extents);
        return {d, data()};
    }

    // Row subscripting.
    Matrix_ref<T, N - 1> operator[](std::size_t n) { return row(n); }

    Matrix_ref<const T, N - 1> operator[](std::size_t n) const
    {
        return row(n);
    }

    // Return a reference to the n'th row of the Matrix_ref.
    Matrix_ref<T, N - 1> row(std::size_t n);
    Matrix_ref<const T, N - 1> row(std::size_t n) const;

    // Return a reference to the n'th column of the Matrix_ref.
    Matrix_ref<T, N - 1> column(std::size_t n);
    Matrix_ref<const T, N - 1> column(std::size_t n) const;

    // Iterators:

    iterator begin() { return {this->desc, ptr}; }
    iterator end() { return {this->desc, ptr, true}; }

    const_iterator begin() const { return {this->desc, ptr}; }
    const_iterator end() const { return {this->desc, ptr, true}; }

    // Mutators:

    void swap(Matrix_ref& m);

    // Apply f(x) for every element x:
    template <typename F>
    Matrix_ref& apply(F f);

    // Apply f(x, mx) for corresponding elements *this and m:
    template <typename M, typename F>
    Enable_if<Matrix_type<M>(), Matrix_ref&> apply(const M& m, F f);

    // Arithmetic operations:

    Matrix_ref& operator=(const T& value); // assignment with scalar

    Matrix_ref& operator+=(const T& value); // scalar addition
    Matrix_ref& operator-=(const T& value); // scalar subtraction
    Matrix_ref& operator*=(const T& value); // scalar multiplication
    Matrix_ref& operator/=(const T& value); // scalar division
    Matrix_ref& operator%=(const T& value); // scalar modulo

    template <typename M>
    Enable_if<Matrix_type<M>(), Matrix_ref&>
    operator+=(const M& m); // matrix addition

    template <typename M>
    Enable_if<Matrix_type<M>(), Matrix_ref&>
    operator-=(const M& m); // matrix subtraction

private:
    T* ptr; // points to the first element of the matrix
};

template <typename T, std::size_t N>
Matrix_ref<T, N>::Matrix_ref(Matrix_ref&& m)
    : Matrix_base<T, N>(m.desc), ptr(m.ptr)
{
}

template <typename T, std::size_t N>
Matrix_ref<T, N>& Matrix_ref<T, N>::operator=(Matrix_ref&& m)
{
    assert(same_extents(this->desc, m.desc));
    std::move(m.begin(), m.end(), begin()); // not exception safe
    return *this;
}

template <typename T, std::size_t N>
Matrix_ref<T, N>::Matrix_ref(const Matrix_ref& m)
    : Matrix_base<T, N>(m.descriptor()), ptr(m.ptr)
{
}

template <typename T, std::size_t N>
template <typename U>
Matrix_ref<T, N>& Matrix_ref<T, N>::operator=(const Matrix_ref<U, N>& m)
{
    static_assert(Convertible<U, T>(),
                  "Matrix_ref assignment: incompatible element types");

    // Just assign values; no resizing is done.
    assert(same_extents(this->desc, m.descriptor()));
    apply(m, [](T& a, const U& b) { a = b; });
    return *this;
}

template <typename T, std::size_t N>
Matrix_ref<T, N>::Matrix_ref(Matrix<value_type, N>& m)
    : Matrix_base<T, N>(m.descriptor()), ptr(m.data())
{
}

template <typename T, std::size_t N>
Matrix_ref<T, N>::Matrix_ref(const Matrix<value_type, N>& m)
    : Matrix_base<T, N>(m.descriptor()), ptr(m.data())
{
}

template <typename T, std::size_t N>
Matrix_ref<T, N>& Matrix_ref<T, N>::operator=(const Matrix<value_type, N>& m)
{
    // Just assign values; no resizing is done.
    assert(same_extents(this->desc, m.descriptor()));
    apply(m, [](T& a, const T& b) { a = b; });
    return *this;
}

template <typename T, std::size_t N>
inline Matrix_ref<T, N - 1> Matrix_ref<T, N>::row(std::size_t n)
{
    assert(n < this->rows());
    auto r = matrix_impl::slice_dim<0>(this->desc, n);
    return {r, ptr};
}

template <typename T, std::size_t N>
inline Matrix_ref<const T, N - 1> Matrix_ref<T, N>::row(std::size_t n) const
{
    assert(n < this->rows());
    auto r = matrix_impl::slice_dim<0>(this->desc, n);
    return {r, ptr};
}

template <typename T, std::size_t N>
inline Matrix_ref<T, N - 1> Matrix_ref<T, N>::column(std::size_t n)
{
    assert(n < this->cols());
    auto c = matrix_impl::slice_dim<1>(this->desc, n);
    return {c, ptr};
}

template <typename T, std::size_t N>
inline Matrix_ref<const T, N - 1> Matrix_ref<T, N>::column(std::size_t n) const
{
    assert(n < this->cols());
    auto c = matrix_impl::slice_dim<1>(this->desc, n);
    return {c, ptr};
}

template <typename T, std::size_t N>
inline void Matrix_ref<T, N>::swap(Matrix_ref& m)
{
    std::swap(this->desc, m.desc);
    std::swap(ptr, m.ptr);
}

template <typename T, std::size_t N>
template <typename F>
Matrix_ref<T, N>& Matrix_ref<T, N>::apply(F f)
{
    for (auto i = begin(); i != end(); ++i) {
        f(*i);
    }
    return *this;
}

template <typename T, std::size_t N>
template <typename M, typename F>
Enable_if<Matrix_type<M>(), Matrix_ref<T, N>&>
Matrix_ref<T, N>::apply(const M& m, F f)
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
inline Matrix_ref<T, N>& Matrix_ref<T, N>::operator=(const T& value)
{
    return apply([&](T& a) { a = value; });
}

template <typename T, std::size_t N>
inline Matrix_ref<T, N>& Matrix_ref<T, N>::operator+=(const T& value)
{
    return apply([&](T& a) { a += value; });
}

template <typename T, std::size_t N>
inline Matrix_ref<T, N>& Matrix_ref<T, N>::operator-=(const T& value)
{
    return apply([&](T& a) { a -= value; });
}

template <typename T, std::size_t N>
inline Matrix_ref<T, N>& Matrix_ref<T, N>::operator*=(const T& value)
{
    return apply([&](T& a) { a *= value; });
}

template <typename T, std::size_t N>
inline Matrix_ref<T, N>& Matrix_ref<T, N>::operator/=(const T& value)
{
    return apply([&](T& a) { a /= value; });
}

template <typename T, std::size_t N>
inline Matrix_ref<T, N>& Matrix_ref<T, N>::operator%=(const T& value)
{
    return apply([&](T& a) { a %= value; });
}

template <typename T, std::size_t N>
template <typename M>
inline Enable_if<Matrix_type<M>(), Matrix_ref<T, N>&>
// clang-format off
Matrix_ref<T, N>::operator+=(const M& m)
// clang-format on
{
#ifdef __clang__ // ugly hack to work around bug in Clang on Mac OS X
    assert(m.order == N);
#else
    static_assert(m.order == N, "+=: mismatched Matrix dimensions");
#endif
    assert(same_extents(this->desc, m.descriptor()));

    return apply(m, [](T& a, const Value_type<M>& b) { a += b; });
}

template <typename T, std::size_t N>
template <typename M>
inline Enable_if<Matrix_type<M>(), Matrix_ref<T, N>&>
// clang-format off
Matrix_ref<T, N>::operator-=(const M& m)
// clang-format on
{
#ifdef __clang__ // ugly hack to work around bug in Clang on Mac OS X
    assert(m.order == N);
#else
    static_assert(m.order == N, "-=: mismatched Matrix dimensions");
#endif
    assert(same_extents(this->desc, m.descriptor()));

    return apply(m, [](T& a, const Value_type<M>& b) { a -= b; });
}

//------------------------------------------------------------------------------

// Specialization:

template <typename T>
class Matrix_ref<T, 0> : public Matrix_base<T, 0> {
public:
    Matrix_ref() = delete;

    Matrix_ref(const Matrix_slice<0>& ms, T* p) : ptr(p + ms.start) {}

    Matrix_ref& operator=(const T& value)
    {
        *ptr = value;
        return *this;
    }

    T& operator()() { return *ptr; }
    const T& operator()() const { return *ptr; }

    operator T&() { return *ptr; }

private:
    T* ptr;
};

} // namespace num

#endif // NUMLIB_MATRIX_MATRIX_REF_H
