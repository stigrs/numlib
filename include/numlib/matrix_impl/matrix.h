// Copyright (c) 2008-2010 Kent State University
// Copyright (c) 2011-2012 Texas A&M University
// Copyright (c) Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#ifndef NUMLIB_MATRIX_MATRIX_H
#define NUMLIB_MATRIX_MATRIX_H

#include <numlib/matrix_impl/matrix_base.h>
#include <initializer_list>
#include <utility>
#include <vector>

namespace Numlib {

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

    // Specify the extents:
    template <typename... Exts>
    explicit Matrix(Exts... exts);

// Construct and assign from Matrix_ref:
#if 0
    template <typename U>
    Matrix(const Matrix_ref<U, N>&);

    template <typename U>
    Matrix& operator=(const Matrix_ref<U, N>&);
#endif
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

    // clang-format off
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
    // clang-format on

    // Iterators:

    iterator begin() { return elems.begin(); }
    iterator end() { return elems.end(); }

    const_iterator begin() const { return elems.begin(); }
    const_iterator end() const { return elems.end(); }

    // Mutators:

    void swap(Matrix& m);

    // Apply f(x) for every element x:
    template <typename F>
    Matrix& apply(F f);

    // Apply f(x, mx) for corresponding elements *this and m:
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
inline Matrix<T, N>::Matrix(Matrix_initializer<T, N> init)
{
    this->desc.start = 0;
    this->desc.extents = Matrix_impl::derive_extents<N>(init);
    Matrix_impl::compute_strides(this->desc);
    elems.reserve(this->desc.size);
    Matrix_impl::insert_flat(init, elems);
    assert(elems.size() == this->desc.size);
}

template <typename T, std::size_t N>
inline Matrix<T, N>& Matrix<T, N>::operator=(Matrix_initializer<T, N> init)
{
    Matrix tmp(init);
    swap(tmp);
    return *this;
}

template <typename T, std::size_t N>
inline void Matrix<T, N>::swap(Matrix& m)
{
    std::swap(this->desc, m.desc);
    elems.swap(m.elems);
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
inline Enable_if<Matrix_type<M>(), Matrix<T, N>&> Matrix<T, N>::
operator+=(const M& m)
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
inline Enable_if<Matrix_type<M>(), Matrix<T, N>&> Matrix<T, N>::
operator-=(const M& m)
{
#ifdef __clang__ // ugly hack to work around bug in Clang on Mac OS X
    assert(m.order == N);
#else
    static_assert(m.order == N, "-=: mismatched Matrix dimensions");
#endif
    assert(same_extents(this->desc, m.descriptor()));

    return apply(m, [](T& a, const Value_type<M>& b) { a -= b; });
}

} // namespace Numlib

#endif // NUMLIB_MATRIX_MATRIX_H
