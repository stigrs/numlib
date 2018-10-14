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

#ifndef SRS_ARRAY_REF3_H
#define SRS_ARRAY_REF3_H

#include <srs/array_impl/functors.h>
#include <srs/types.h>
#include <array>
#include <gsl/gsl>


namespace srs {

//
// Three-dimensional dense array reference.
//
template <class T>
class Array_ref<T, 3> {
public:
    static constexpr int rank = 3;

    typedef T value_type;
    typedef Int_t size_type;

    // Constructors:
    Array_ref(size_type n1,
              size_type n2,
              size_type n3,
              size_type str1,
              size_type str2,
              T* ptr)
        : elems(ptr), extents{n1, n2, n3}, strides{str1, str2}
    {
    }

    // Assignment:

    Array_ref& operator=(const Array<T, 3>& a);

    // Element access:

    T& at(size_type i, size_type j, size_type k);
    const T& at(size_type i, size_type j, size_type k) const;

    T& operator()(size_type i, size_type j, size_type k);
    const T& operator()(size_type i, size_type j, size_type k) const;

    Array_ref<T, 2> operator[](size_type i) { return depth(i); }
    Array_ref<const T, 2> operator[](size_type i) const { return depth(i); }

    // Slicing:

    Array_ref<T, 2> depth(size_type i);
    Array_ref<const T, 2> depth(size_type i) const;

    // Capacity:

    size_type rows() const { return extents[0]; }
    size_type cols() const { return extents[1]; }
    size_type depths() const { return extents[2]; }
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
    std::array<size_type, 3> extents;
    std::array<size_type, 2> strides;
};

template <class T>
Array_ref<T, 3>& Array_ref<T, 3>::operator=(const Array<T, 3>& a)
{
    Expects(rows() == a.rows());
    Expects(cols() == a.cols());
    Expects(depths() == a.depths());
    for (size_type k = 0; k < extents[2]; ++k) {
        for (size_type j = 0; j < extents[1]; ++j) {
            for (size_type i = 0; i < extents[0]; ++i) {
                (*this)(i, j, k) = a(i, j, k);
            }
        }
    }
    return *this;
}

template <class T>
inline T& Array_ref<T, 3>::at(size_type i, size_type j, size_type k)
{
    Expects(i >= 0 && i < extents[0]);
    Expects(j >= 0 && j < extents[1]);
    Expects(k >= 0 && k < extents[2]);
    return elems[i + j * strides[0] + k * strides[1]];
}

template <class T>
inline const T& Array_ref<T, 3>::at(size_type i, size_type j, size_type k) const
{
    Expects(i >= 0 && i < extents[0]);
    Expects(j >= 0 && j < extents[1]);
    Expects(k >= 0 && k < extents[2]);
    return elems[i + j * strides[0] + k * strides[1]];
}

template <class T>
inline T& Array_ref<T, 3>::operator()(size_type i, size_type j, size_type k)
{
#ifdef NDEBUG
    return elems[i + j * strides[0] + k * strides[1]];
#else
    return at(i, j, k);
#endif
}

template <class T>
inline const T& Array_ref<T, 3>::operator()(size_type i,
                                            size_type j,
                                            size_type k) const
{
#ifdef NDEBUG
    return elems[i + j * strides[0] + k * strides[1]];
#else
    return at(i, j, k);
#endif
}

template <class T>
inline Array_ref<T, 2> Array_ref<T, 3>::depth(size_type i)
{
    Expects(i >= 0 && i < extents[2]);
    return Array_ref<T, 2>(
        extents[0], extents[1], strides[0], data() + i * strides[1]);
}

template <class T>
inline Array_ref<const T, 2> Array_ref<T, 3>::depth(size_type i) const
{
    Expects(i >= 0 && i < extents[2]);
    return Array_ref<const T, 2>(
        extents[0], extents[1], strides[0], data() + i * strides[1]);
}

template <class T>
template <class F>
inline Array_ref<T, 3>& Array_ref<T, 3>::apply(F f)
{
    for (size_type k = 0; k < extents[2]; ++k) {
        for (size_type j = 0; j < extents[1]; ++j) {
            for (size_type i = 0; i < extents[0]; ++i) {
                f((*this)(i, j, k));
            }
        }
    }
    return *this;
}

template <class T>
template <class F>
inline Array_ref<T, 3>& Array_ref<T, 3>::apply(F f, const T& value)
{
    for (size_type k = 0; k < extents[2]; ++k) {
        for (size_type j = 0; j < extents[1]; ++j) {
            for (size_type i = 0; i < extents[0]; ++i) {
                f((*this)(i, j, k), value);
            }
        }
    }
    return *this;
}

template <class T>
inline Array_ref<T, 3>& Array_ref<T, 3>::operator=(const T& value)
{
    return apply(Assign<T>(), value);
}

template <class T>
inline Array_ref<T, 3>& Array_ref<T, 3>::operator*=(const T& value)
{
    return apply(Mul_assign<T>(), value);
}

template <class T>
inline Array_ref<T, 3>& Array_ref<T, 3>::operator/=(const T& value)
{
    return apply(Div_assign<T>(), value);
}

template <class T>
inline Array_ref<T, 3>& Array_ref<T, 3>::operator%=(const T& value)
{
    return apply(Mod_assign<T>(), value);
}

template <class T>
inline Array_ref<T, 3>& Array_ref<T, 3>::operator+=(const T& value)
{
    return apply(Add_assign<T>(), value);
}

template <class T>
inline Array_ref<T, 3>& Array_ref<T, 3>::operator-=(const T& value)
{
    return apply(Minus_assign<T>(), value);
}

template <class T>
inline Array_ref<T, 3>& Array_ref<T, 3>::operator&=(const T& value)
{
    return apply(And_assign<T>(), value);
}

template <class T>
inline Array_ref<T, 3>& Array_ref<T, 3>::operator|=(const T& value)
{
    return apply(Or_assign<T>(), value);
}

template <class T>
inline Array_ref<T, 3>& Array_ref<T, 3>::operator^=(const T& value)
{
    return apply(Xor_assign<T>(), value);
}

template <class T>
inline Array_ref<T, 3>& Array_ref<T, 3>::operator!()
{
    return apply(Not<T>());
}

template <class T>
inline Array_ref<T, 3>& Array_ref<T, 3>::operator-()
{
    return apply(Unary_minus<T>());
}

template <class T>
inline Array_ref<T, 3>& Array_ref<T, 3>::operator~()
{
    return apply(Complement<T>());
}

template <class T>
Array_ref<T, 3>& Array_ref<T, 3>::operator+=(const Array_ref<T, 3>& a)
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
Array_ref<T, 3>& Array_ref<T, 3>::operator-=(const Array_ref<T, 3>& a)
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

#endif  // SRS_ARRAY_REF3_H
