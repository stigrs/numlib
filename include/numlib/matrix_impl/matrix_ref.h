// Copyright (c) 2008-2010 Kent State University
// Copyright (c) 2011-2012 Texas A&M University
// Copyright (c) Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#ifndef NUMLIB_MATRIX_MATRIX_REF_H
#define NUMLIB_MATRIX_MATRIX_REF_H

#include <numlib/matrix_impl/matrix_base.h>

namespace Numlib {

//------------------------------------------------------------------------------

template <typename T, std::size_t N>
class Matrix_ref : public Matrix_base<T, N> {
public:
    using size_type = typename Matrix_base<T, N>::size_type;
    using value_type = Remove_const<T>;

    Matrix_ref() = default;

    // Move construction and assignment:
    Matrix_ref(Matrix_ref&&) = default;
    Matrix_ref& operator=(Matrix_ref&&) = default;

    // Copy construction and assignment:
    Matrix_ref(const Matrix_ref&) = default;
    Matrix_ref& operator=(const Matrix_ref&) = default;

    // Specify the extents:
    template <typename... Exts>
    explicit Matrix(Exts... exts);

    // Construct and assign from Matrix_ref:
    template <typename U>
    Matrix(const Matrix_ref<U, N>&);
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

} // namespace Numlib

#endif // NUMLIB_MATRIX_MATRIX_REF_H
