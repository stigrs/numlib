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
    Matrix_ref(Matrix_ref&&);
    Matrix_ref& operator=(Matrix_ref&&);

    // Copy construction and assignment:
    Matrix_ref(const Matrix_ref&);
    Matrix_ref& operator=(const Matrix_ref&);

    // Slice initialization:
    //
    // Initialize the Matrix_ref over the Matrix_slice, starting at the
    // element pointed to by p.
    Matrix_ref(const Matrix_slice<N>& ms, T* p);

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

    ~Matrix() = default;

    // "Flat" element access:
    T* data() { return ptr; }
    const T* data() const { return ptr; }

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

    // Mutators:

    void swap(Matrix_ref& m);

private:
    T* ptr; // points to the first element of the matrix
};

template <typename T, std::size_t N>
inline void Matrix_ref<T, N>::swap(Matrix_ref& m)
{
    std::swap(this->desc, m.desc);
    std::swap(ptr, m.ptr);
}

} // namespace Numlib

#endif // NUMLIB_MATRIX_MATRIX_REF_H
