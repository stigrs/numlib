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

#ifndef SRS_ARRAY_REF1_H
#define SRS_ARRAY_REF1_H

#include <srs/array.h>
#include <srs/array_impl/functors.h>
#include <srs/array_impl/slice_iter.h>
#include <srs/types.h>
#include <array>
#include <gsl/gsl>


namespace srs {

//
// One-dimensional dense array reference.
//
template <class T>
class Array_ref<T, 1> {
public:
    static constexpr int rank = 1;

    typedef T value_type;
    typedef Int_t size_type;
    typedef Slice_iter<T> iterator;
    typedef Slice_iter<const T> const_iterator;

    // Constructors:

    Array_ref(size_type n, size_type str, T* ptr)
        : elems(ptr), extents{n}, stride(str)
    {
    }

    // Assignment:
    Array_ref& operator=(const Array<T, 1>& a);

    // Element access:

    T& at(size_type i);
    const T& at(size_type i) const;

    T& operator()(size_type i);
    const T& operator()(size_type i) const;

    T& operator[](size_type i);
    const T& operator[](size_type i) const;

    // Iterators:

    iterator begin();
    iterator end();

    const_iterator begin() const;
    const_iterator end() const;

    // Capacity:

    size_type size() const { return extents[0]; }

    // Access underlying array:

    T* data() { return elems; }
    const T* data() const { return elems; }

    // Modifiers:

    void swap(const Array_ref& a);

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
    std::array<size_type, 1> extents;
    size_type stride;
};

template <class T>
inline Array_ref<T, 1>& Array_ref<T, 1>::operator=(const Array<T, 1>& a)
{
    Expects(size() == a.size());
    for (size_type i = 0; i < extents[0]; ++i) {
        (*this)(i) = a(i);
    }
    return *this;
}

template <class T>
inline T& Array_ref<T, 1>::at(size_type i)
{
    Expects(i >= 0 && i < extents[0]);
    return elems[i * stride];
}

template <class T>
inline const T& Array_ref<T, 1>::at(size_type i) const
{
    Expects(i >= 0 && i < extents[0]);
    return elems[i * stride];
}

template <class T>
inline T& Array_ref<T, 1>::operator()(size_type i)
{
#ifdef NDEBUG
    return elems[i * stride];
#else
    return at(i);
#endif
}

template <class T>
inline const T& Array_ref<T, 1>::operator()(size_type i) const
{
#ifdef NDEBUG
    return elems[i * stride];
#else
    return at(i);
#endif
}

template <class T>
inline T& Array_ref<T, 1>::operator[](size_type i)
{
#ifdef NDEBUG
    return elems[i * stride];
#else
    return at(i);
#endif
}

template <class T>
inline const T& Array_ref<T, 1>::operator[](size_type i) const
{
#ifdef NDEBUG
    return elems[i * stride];
#else
    return at(i);
#endif
}

template <class T>
Slice_iter<T> Array_ref<T, 1>::begin()
{
    return Slice_iter<T>(elems, srs::Slice(0, size(), stride));
}

template <class T>
Slice_iter<const T> Array_ref<T, 1>::begin() const
{
    return Slice_iter<const T>(elems, srs::Slice(0, size(), stride));
}

template <class T>
Slice_iter<T> Array_ref<T, 1>::end()
{
    return Slice_iter<T>(elems, srs::Slice(size(), 1, stride));
}

template <class T>
Slice_iter<const T> Array_ref<T, 1>::end() const
{
    return Slice_iter<const T>(elems, srs::Slice(size(), 1, stride));
}

template <class T>
inline void Array_ref<T, 1>::swap(const Array_ref<T, 1>& a)
{
    Expects(extents = a.extents);
    Expects(stride = a.stride);
    for (size_type i = 0; i < extents[0]; ++i) {
        std::swap(elems[i], a(i));
    }
}

template <class T>
template <class F>
inline Array_ref<T, 1>& Array_ref<T, 1>::apply(F f)
{
    for (size_type i = 0; i < extents[0]; ++i) {
        f((*this)(i));
    }
    return *this;
}

template <class T>
template <class F>
inline Array_ref<T, 1>& Array_ref<T, 1>::apply(F f, const T& value)
{
    for (size_type i = 0; i < extents[0]; ++i) {
        f((*this)(i), value);
    }
    return *this;
}

template <class T>
inline Array_ref<T, 1>& Array_ref<T, 1>::operator=(const T& value)
{
    return apply(Assign<T>(), value);
}

template <class T>
inline Array_ref<T, 1>& Array_ref<T, 1>::operator*=(const T& value)
{
    return apply(Mul_assign<T>(), value);
}

template <class T>
inline Array_ref<T, 1>& Array_ref<T, 1>::operator/=(const T& value)
{
    return apply(Div_assign<T>(), value);
}

template <class T>
inline Array_ref<T, 1>& Array_ref<T, 1>::operator%=(const T& value)
{
    return apply(Mod_assign<T>(), value);
}

template <class T>
inline Array_ref<T, 1>& Array_ref<T, 1>::operator+=(const T& value)
{
    return apply(Add_assign<T>(), value);
}

template <class T>
inline Array_ref<T, 1>& Array_ref<T, 1>::operator-=(const T& value)
{
    return apply(Minus_assign<T>(), value);
}

template <class T>
inline Array_ref<T, 1>& Array_ref<T, 1>::operator&=(const T& value)
{
    return apply(And_assign<T>(), value);
}

template <class T>
inline Array_ref<T, 1>& Array_ref<T, 1>::operator|=(const T& value)
{
    return apply(Or_assign<T>(), value);
}

template <class T>
inline Array_ref<T, 1>& Array_ref<T, 1>::operator^=(const T& value)
{
    return apply(Xor_assign<T>(), value);
}

template <class T>
inline Array_ref<T, 1>& Array_ref<T, 1>::operator!()
{
    return apply(Not<T>());
}

template <class T>
inline Array_ref<T, 1>& Array_ref<T, 1>::operator-()
{
    return apply(Unary_minus<T>());
}

template <class T>
inline Array_ref<T, 1>& Array_ref<T, 1>::operator~()
{
    return apply(Complement<T>());
}

template <class T>
Array_ref<T, 1>& Array_ref<T, 1>::operator+=(const Array_ref<T, 1>& a)
{
    Expects(extents == a.extents);
    for (size_type i = 0; i < a.size(); ++i) {
        (*this)(i) += a(i);
    }
    return *this;
}

template <class T>
Array_ref<T, 1>& Array_ref<T, 1>::operator-=(const Array_ref<T, 1>& a)
{
    Expects(extents == a.extents);
    for (size_type i = 0; i < a.size(); ++i) {
        (*this)(i) -= a(i);
    }
    return *this;
}

}  // namespace srs

#endif  // SRS_ARRAY_REF1_H
