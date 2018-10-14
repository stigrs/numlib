////////////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2017 Stig Rune Sellevag. All rights reserved.
//
// This code is licensed under the MIT License (MIT).
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.
//
////////////////////////////////////////////////////////////////////////////////

#ifndef SRS_ARRAY_REF2_H
#define SRS_ARRAY_REF2_H

#include <srs/array_impl/functors.h>
#include <srs/types.h>
#include <array>
#include <gsl/gsl>


namespace srs {

//
// Two-dimensional dense array reference.
//
template <class T>
class Array_ref<T, 2> {
public:
    static constexpr int rank = 2;

    typedef T value_type;
    typedef Int_t size_type;

    // Constructors:

    Array_ref(size_type nrows, size_type ncols, size_type str, T* ptr)
        : elems(ptr), extents{nrows, ncols}, stride(str)
    {
    }

    // Assignment:

    Array_ref& operator=(const Array<T, 2>& a);

    // Element access:

    T& at(size_type i, size_type j);
    const T& at(size_type i, size_type j) const;

    T& operator()(size_type i, size_type j);
    const T& operator()(size_type i, size_type j) const;

    Array_ref<T, 1> operator[](size_type i) { return row(i); }
    Array<T, 1> operator[](size_type i) const { return row(i); }

    // Slicing:

    Array_ref<T, 1> row(size_type i);
    Array_ref<const T, 1> row(size_type i) const;

    Array_ref<T, 1> column(size_type i);
    Array_ref<const T, 1> column(size_type i) const;

    Array_ref<T, 1> diag();
    Array_ref<const T, 1> diag() const;

	// Flatten matrix to one-dimensional:

	Array_ref<T, 1> flatten();
	Array_ref<const T, 1> flatten() const;


    // Capacity:

	size_type size() const { return extents[0] * extents[1]; }
    size_type rows() const { return extents[0]; }
    size_type cols() const { return extents[1]; }
    size_type extent(size_type dim) const
    {
        Expects(dim >= 0 && dim < rank);
        return extents[dim];
    }

    // Access underlying array:

    T* data() { return elems; }
    const T* data() const { return elems; }

    // Element-wise operations:

    template <class F>
    Array_ref& apply(F f);

    template <class F>
    Array_ref& apply(F f, const T& value);

    Array_ref& operator=(const T& value);

    Array_ref& operator*=(const T& value);
    Array_ref& operator/=(const T& value);
    Array_ref& operator%=(const T& value);
    Array_ref& operator+=(const T& value);
    Array_ref& operator-=(const T& value);

    Array_ref& operator&=(const T& value);
    Array_ref& operator|=(const T& value);
    Array_ref& operator^=(const T& value);

    Array_ref& operator!();
    Array_ref& operator-();
    Array_ref& operator~();

    Array_ref& operator+=(const Array_ref& a);
    Array_ref& operator-=(const Array_ref& a);

private:
    T* elems;
    std::array<size_type, 2> extents;
    size_type stride;
};

template <class T>
Array_ref<T, 2>& Array_ref<T, 2>::operator=(const Array<T, 2>& a)
{
    Expects(rows() == a.rows());
    Expects(cols() == a.cols());
    for (size_type j = 0; j < extents[1]; ++j) {
        for (size_type i = 0; i < extents[0]; ++i) {
            (*this)(i, j) = a(i, j);
        }
    }
    return *this;
}

template <class T>
inline T& Array_ref<T, 2>::at(size_type i, size_type j)
{
    Expects(i >= 0 && i < extents[0] && j >= 0 && j < extents[1]);
    return elems[i + j * stride];
}

template <class T>
inline const T& Array_ref<T, 2>::at(size_type i, size_type j) const
{
    Expects(i >= 0 && i < extents[0] && j >= 0 && j < extents[1]);
    return elems[i + j * stride];
}

template <class T>
inline T& Array_ref<T, 2>::operator()(size_type i, size_type j)
{
#ifdef NDEBUG
    return elems[i + j * stride];
#else
    return at(i, j);
#endif
}

template <class T>
inline const T& Array_ref<T, 2>::operator()(size_type i, size_type j) const
{
#ifdef NDEBUG
    return elems[i + j * stride];
#else
    return at(i, j);
#endif
}

template <class T>
inline Array_ref<T, 1> Array_ref<T, 2>::row(size_type i)
{
    Expects(i >= 0 && i < extents[0]);
    return Array_ref<T, 1>(extents[1], stride, data() + i);
}

template <class T>
inline Array_ref<const T, 1> Array_ref<T, 2>::row(size_type i) const
{
    Expects(i >= 0 && i < extents[0]);
    return Array_ref<const T, 1>(extents[1], stride, data() + i);
}

template <class T>
inline Array_ref<T, 1> Array_ref<T, 2>::column(size_type i)
{
    Expects(i >= 0 && i < extents[1]);
    return Array_ref<T, 1>(extents[0], 1, data() + i * stride);
}

template <class T>
inline Array_ref<const T, 1> Array_ref<T, 2>::column(size_type i) const
{
    Expects(i >= 0 && i < extents[1]);
    return Array_ref<const T, 1>(extents[0], 1, data() + i * stride);
}

template <class T>
inline Array_ref<T, 1> Array_ref<T, 2>::diag()
{
    Expects(extents[0] == extents[1]);
    return Array_ref<T, 1>(extents[0], stride + 1, data());
}

template <class T>
inline Array_ref<const T, 1> Array_ref<T, 2>::diag() const
{
    Expects(extents[0] == extents[1]);
    return Array_ref<const T, 1>(extents[0], stride + 1, data());
}

template <class T>
inline Array_ref<T, 1> Array_ref<T, 2>::flatten()
{
	return Array_ref<T, 1>(size(), 1, data());
}

template <class T>
inline Array_ref<const T, 1> Array_ref<T, 2>::flatten() const
{ 
	return Array_ref<const T, 1>(size(), 1, data());
}

template <class T>
template <class F>
inline Array_ref<T, 2>& Array_ref<T, 2>::apply(F f)
{
    for (size_type j = 0; j < extents[1]; ++j) {
        for (size_type i = 0; i < extents[0]; ++i) {
            f((*this)(i, j));
        }
    }
    return *this;
}

template <class T>
template <class F>
inline Array_ref<T, 2>& Array_ref<T, 2>::apply(F f, const T& value)
{
    for (size_type j = 0; j < extents[1]; ++j) {
        for (size_type i = 0; i < extents[0]; ++i) {
            f((*this)(i, j), value);
        }
    }
    return *this;
}

template <class T>
inline Array_ref<T, 2>& Array_ref<T, 2>::operator=(const T& value)
{
    return apply(Assign<T>(), value);
}

template <class T>
inline Array_ref<T, 2>& Array_ref<T, 2>::operator*=(const T& value)
{
    return apply(Mul_assign<T>(), value);
}

template <class T>
inline Array_ref<T, 2>& Array_ref<T, 2>::operator/=(const T& value)
{
    return apply(Div_assign<T>(), value);
}

template <class T>
inline Array_ref<T, 2>& Array_ref<T, 2>::operator%=(const T& value)
{
    return apply(Mod_assign<T>(), value);
}

template <class T>
inline Array_ref<T, 2>& Array_ref<T, 2>::operator+=(const T& value)
{
    return apply(Add_assign<T>(), value);
}

template <class T>
inline Array_ref<T, 2>& Array_ref<T, 2>::operator-=(const T& value)
{
    return apply(Minus_assign<T>(), value);
}

template <class T>
inline Array_ref<T, 2>& Array_ref<T, 2>::operator&=(const T& value)
{
    return apply(And_assign<T>(), value);
}

template <class T>
inline Array_ref<T, 2>& Array_ref<T, 2>::operator|=(const T& value)
{
    return apply(Or_assign<T>(), value);
}

template <class T>
inline Array_ref<T, 2>& Array_ref<T, 2>::operator^=(const T& value)
{
    return apply(Xor_assign<T>(), value);
}

template <class T>
inline Array_ref<T, 2>& Array_ref<T, 2>::operator!()
{
    return apply(Not<T>());
}

template <class T>
inline Array_ref<T, 2>& Array_ref<T, 2>::operator-()
{
    return apply(Unary_minus<T>());
}

template <class T>
inline Array_ref<T, 2>& Array_ref<T, 2>::operator~()
{
    return apply(Complement<T>());
}

template <class T>
Array_ref<T, 2>& Array_ref<T, 2>::operator+=(const Array_ref<T, 2>& a)
{
    Expects(extents == a.extents);
    for (size_type j = 0; j < a.cols(); ++j) {
        for (size_type i = 0; i < a.rows(); ++i) {
            (*this)(i, j) += a(i, j);
        }
    }
    return *this;
}

template <class T>
Array_ref<T, 2>& Array_ref<T, 2>::operator-=(const Array_ref<T, 2>& a)
{
    Expects(extents == a.extents);
    for (size_type j = 0; j < a.cols(); ++j) {
        for (size_type i = 0; i < a.rows(); ++i) {
            (*this)(i, j) -= a(i, j);
        }
    }
    return *this;
}

}  // namespace srs

#endif  // SRS_ARRAY_REF2_H
