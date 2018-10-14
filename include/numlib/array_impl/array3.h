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

#ifndef SRS_ARRAY3_H
#define SRS_ARRAY3_H

#include <srs/array_impl/array_ref.h>
#include <srs/array_impl/functors.h>
#include <srs/types.h>
#include <algorithm>
#include <array>
#include <gsl/gsl>
#include <initializer_list>
#include <vector>


namespace srs {

//
// Three-dimensional dense array (cube) class.
//
template <class T>
class Array<T, 3> {
public:
    static constexpr int rank = 3;

    typedef T value_type;
    typedef Int_t size_type;
    typedef typename std::vector<T>::iterator iterator;
    typedef typename std::vector<T>::const_iterator const_iterator;
    typedef typename std::
        initializer_list<std::initializer_list<std::initializer_list<T>>>
            initializer_list_3d;

    // Constructors:

    Array() : elems(), extents{0, 0, 0}, strides{0, 0} {}

    Array(size_type n1, size_type n2, size_type n3)
        : elems(n1 * n2 * n3), extents{n1, n2, n3}, strides{n1, n1 * n2}
    {
    }

    Array(size_type n1, size_type n2, size_type n3, const T& value)
        : elems(n1 * n2 * n3, value), extents{n1, n2, n3}, strides{n1, n1 * n2}
    {
    }

    Array(size_type n1, size_type n2, size_type n3, T* ptr);

    template <Int_t n1, Int_t n2, Int_t n3>
    Array(const T (&a)[n1][n2][n3]);

    Array(initializer_list_3d ilist) { assign(ilist); }

    // Copy elements referenced by array slice.
    template <class U>
    Array(const Array_ref<U, 3>& a);

    // T f(const T&) would be a typical type for f.
    template <class F>
    Array(const Array& a, F f);

    // T f(const T&, const Arg&) would be a typical type for f.
    template <class F, class Arg>
    Array(const Array& a, F f, const Arg& value);

    // Assignments:

    template <class U>
    Array& operator=(const Array_ref<U, 3>& a);
    Array& operator=(initializer_list_3d ilist);

    // Element access:

    T& at(size_type i, size_type j, size_type k);
    const T& at(size_type i, size_type j, size_type k) const;

    T& operator()(size_type i, size_type j, size_type k);
    const T& operator()(size_type i, size_type j, size_type k) const;

    Array_ref<T, 2> operator[](size_type i) { return depth(i); }
    Array_ref<const T, 2> operator[](size_type i) const { return depth(i); }

    // Iterators:

    iterator begin() { return elems.begin(); }
    iterator end() { return elems.end(); }

    const_iterator begin() const { return elems.begin(); }
    const_iterator end() const { return elems.end(); }

    // Slicing:

    Array_ref<T, 2> depth(size_type i);
    Array_ref<const T, 2> depth(size_type i) const;

    Array_ref<T, 3> slice(size_type ifirst,
                          size_type ilast,
                          size_type jfirst,
                          size_type jlast,
                          size_type kfirst,
                          size_type klast);
    Array_ref<const T, 3> slice(size_type ifirst,
                                size_type ilast,
                                size_type jfirst,
                                size_type jlast,
                                size_type kfirst,
                                size_type klast) const;

    // Capacity:

    bool empty() const { return elems.empty(); }

    size_type size() const { return elems.size(); }
    size_type max_size() const { return elems.max_size(); }
    size_type capacity() const { return elems.capacity(); }
    size_type rows() const { return extents[0]; }
    size_type cols() const { return extents[1]; }
    size_type depths() const { return extents[2]; }
    size_type dim1() const { return extents[0]; }
    size_type dim2() const { return extents[1]; }
    size_type dim3() const { return extents[2]; }
    size_type extent(size_type dim) const
    {
        Expects(dim >= 0 && dim < rank);
        return extents[dim];
    }

    // Modifiers:

    void clear();

    void swap(Array& a);

    void resize(size_type n1, size_type n2, size_type n3);
    void resize(size_type n1, size_type n2, size_type n3, const T& value);

    // Access underlying array:

    T* data() { return elems.data(); }
    const T* data() const { return elems.data(); }

    // Element-wise operations:

    template <class F>
    Array& apply(F f);

    template <class F>
    Array& apply(F f, const T& value);

    Array& operator=(const T& value);

    Array& operator*=(const T& value);
    Array& operator/=(const T& value);
    Array& operator%=(const T& value);
    Array& operator+=(const T& value);
    Array& operator-=(const T& value);

    Array& operator&=(const T& value);
    Array& operator|=(const T& value);
    Array& operator^=(const T& value);

    Array& operator!();
    Array& operator-();
    Array& operator~();

    Array& operator+=(const Array& a);
    Array& operator-=(const Array& a);

private:
    std::vector<T> elems;  // storage
    std::array<size_type, 3> extents;
    std::array<size_type, 2> strides;

    // Compute index.
    Int_t index(size_type i, size_type j, size_type k) const;

    // Helper function for assigning initializer list.
    void assign(initializer_list_3d ilist);
};

template <class T>
Array<T, 3>::Array(size_type n1, size_type n2, size_type n3, T* ptr)
    : elems(n1 * n2 * n3), extents{n1, n2, n3}, strides{n1, n1 * n2}
{
    for (size_type i = 0; i < size(); ++i) {
        elems[i] = ptr[i];
    }
}

template <class T>
template <Int_t n1, Int_t n2, Int_t n3>
Array<T, 3>::Array(const T (&a)[n1][n2][n3])
    : elems(n1 * n2 * n3), extents{n1, n2, n3}, strides{n1, n1 * n2}
{
    for (size_type i = 0; i < extents[0]; ++i) {
        for (size_type j = 0; j < extents[1]; ++j) {
            for (size_type k = 0; k < extents[2]; ++k) {
                (*this)(i, j, k) = a[i][j][k];
            }
        }
    }
}

template <class T>
template <class U>
Array<T, 3>::Array(const Array_ref<U, 3>& a)
    : elems(a.rows() * a.cols() * a.depths()),
      extents{a.rows(), a.cols(), a.depths()},
      strides{a.rows(), a.rows() * a.cols()}
{
    for (size_type k = 0; k < a.depths(); ++k) {
        for (size_type j = 0; j < a.cols(); ++j) {
            for (size_type i = 0; i < a.rows(); ++i) {
                (*this)(i, j, k) = a(i, j, k);
            }
        }
    }
}

template <class T>
template <class F>
Array<T, 3>::Array(const Array& a, F f)
{
    Expects(size() == a.size());
    for (size_type i = 0; i < size(); ++i) {
        elems[i] = f(a.data()[i]);
    }
}

template <class T>
template <class F, class Arg>
Array<T, 3>::Array(const Array& a, F f, const Arg& value)
{
    Expects(size() == a.size());
    for (size_type i = 0; i < size(); ++i) {
        elems[i] = f(a.data()[i], value);
    }
}

template <class T>
template <class U>
Array<T, 3>& Array<T, 3>::operator=(const Array_ref<U, 3>& a)
{
    resize(a.rows() * a.cols() * a.depths());
    extents = {a.rows(), a.cols(), a.depths()};
    strides = {a.rows(), a.rows() * a.cols()};

    for (size_type k = 0; k < a.depths(); ++k) {
        for (size_type j = 0; j < a.cols(); ++j) {
            for (size_type i = 0; i < a.rows(); ++i) {
                (*this)(i, j, k) = a(i, j, k);
            }
        }
    }
    return *this;
}

template <class T>
inline Array<T, 3>& Array<T, 3>::operator=(initializer_list_3d ilist)
{
    assign(ilist);
    return *this;
}

template <class T>
inline T& Array<T, 3>::at(size_type i, size_type j, size_type k)
{
    Expects(i >= 0 && i < extents[0]);
    Expects(j >= 0 && j < extents[1]);
    Expects(k >= 0 && k < extents[2]);
    return elems[index(i, j, k)];
}

template <class T>
inline const T& Array<T, 3>::at(size_type i, size_type j, size_type k) const
{
    Expects(i >= 0 && i < extents[0]);
    Expects(j >= 0 && j < extents[1]);
    Expects(k >= 0 && k < extents[2]);
    return elems[index(i, j, k)];
}

template <class T>
inline T& Array<T, 3>::operator()(size_type i, size_type j, size_type k)
{
#ifdef NDEBUG
    return elems[index(i, j, k)];
#else
    return at(i, j, k);
#endif
}

template <class T>
inline const T& Array<T, 3>::operator()(size_type i,
                                        size_type j,
                                        size_type k) const
{
#ifdef NDEBUG
    return elems[index(i, j, k)];
#else
    return at(i, j, k);
#endif
}

template <class T>
inline Array_ref<T, 2> Array<T, 3>::depth(size_type i)
{
    Expects(i >= 0 && i < extents[2]);
    return Array_ref<T, 2>(
        extents[0], extents[1], strides[0], data() + i * strides[1]);
}

template <class T>
inline Array_ref<const T, 2> Array<T, 3>::depth(size_type i) const
{
    Expects(i >= 0 && i < extents[2]);
    return Array_ref<const T, 2>(
        extents[0], extents[1], strides[0], data() + i * strides[1]);
}

template <class T>
Array_ref<T, 3> Array<T, 3>::slice(size_type ifirst,
                                   size_type ilast,
                                   size_type jfirst,
                                   size_type jlast,
                                   size_type kfirst,
                                   size_type klast)
{
    Expects(ifirst >= 0 && ifirst < ilast && ilast < extents[0]);
    Expects(jfirst >= 0 && jfirst < jlast && jlast < extents[1]);
    Expects(kfirst >= 0 && kfirst < klast && klast < extents[2]);
    size_type n1 = ilast - ifirst + 1;
    size_type n2 = jlast - jfirst + 1;
    size_type n3 = klast - kfirst + 1;
    // clang-format off
    return Array_ref<T, 3>(n1, n2, n3, strides[0], strides[1], 
        data() + ifirst + jfirst * strides[0] + kfirst * strides[1]);
    // clang-format on
}

template <class T>
Array_ref<const T, 3> Array<T, 3>::slice(size_type ifirst,
                                         size_type ilast,
                                         size_type jfirst,
                                         size_type jlast,
                                         size_type kfirst,
                                         size_type klast) const
{
    Expects(ifirst >= 0 && ifirst < ilast && ilast < extents[0]);
    Expects(jfirst >= 0 && jfirst < jlast && jlast < extents[1]);
    Expects(kfirst >= 0 && kfirst < klast && klast < extents[2]);
    size_type n1 = ilast - ifirst + 1;
    size_type n2 = jlast - jfirst + 1;
    size_type n3 = klast - kfirst + 1;
    // clang-format off
    return Array_ref<const T, 3>(n1, n2, n3, strides[0], strides[1], 
        data() + ifirst + jfirst * strides[0] + kfirst * strides[1]);
    // clang-format on
}

template <class T>
inline void Array<T, 3>::clear()
{
    elems.clear();
    extents = {0, 0, 0};
    strides = {0, 0};
}

template <class T>
inline void Array<T, 3>::swap(Array<T, 3>& a)
{
    elems.swap(a.elems);
    std::swap(extents, a.extents);
    std::swap(strides, a.strides);
}

template <class T>
inline void Array<T, 3>::resize(size_type n1, size_type n2, size_type n3)
{
    elems.resize(n1 * n2 * n3);
    extents = {n1, n2, n3};
    strides = {n1, n1 * n2};
}

template <class T>
inline void Array<T, 3>::resize(size_type n1,
                                size_type n2,
                                size_type n3,
                                const T& value)
{
    elems.resize(n1 * n2 * n3, value);
    extents = {n1, n2, n3};
    strides = {n1, n1 * n2};
}

template <class T>
template <class F>
inline Array<T, 3>& Array<T, 3>::apply(F f)
{
    for (auto& v : elems) {
        f(v);
    }
    return *this;
}

template <class T>
template <class F>
inline Array<T, 3>& Array<T, 3>::apply(F f, const T& value)
{
    for (auto& v : elems) {
        f(v, value);
    }
    return *this;
}

template <class T>
inline Array<T, 3>& Array<T, 3>::operator=(const T& value)
{
    apply(Assign<T>(), value);
    return *this;
}

template <class T>
inline Array<T, 3>& Array<T, 3>::operator*=(const T& value)
{
    apply(Mul_assign<T>(), value);
    return *this;
}

template <class T>
inline Array<T, 3>& Array<T, 3>::operator/=(const T& value)
{
    apply(Div_assign<T>(), value);
    return *this;
}

template <class T>
inline Array<T, 3>& Array<T, 3>::operator%=(const T& value)
{
    apply(Mod_assign<T>(), value);
    return *this;
}

template <class T>
inline Array<T, 3>& Array<T, 3>::operator+=(const T& value)
{
    apply(Add_assign<T>(), value);
    return *this;
}

template <class T>
inline Array<T, 3>& Array<T, 3>::operator-=(const T& value)
{
    apply(Minus_assign<T>(), value);
    return *this;
}

template <class T>
inline Array<T, 3>& Array<T, 3>::operator&=(const T& value)
{
    apply(And_assign<T>(), value);
    return *this;
}

template <class T>
inline Array<T, 3>& Array<T, 3>::operator|=(const T& value)
{
    apply(Or_assign<T>(), value);
    return *this;
}

template <class T>
inline Array<T, 3>& Array<T, 3>::operator^=(const T& value)
{
    apply(Xor_assign<T>(), value);
    return *this;
}

template <class T>
inline Array<T, 3>& Array<T, 3>::operator!()
{
    apply(Not<T>());
    return *this;
}

template <class T>
inline Array<T, 3>& Array<T, 3>::operator-()
{
    apply(Unary_minus<T>());
    return *this;
}

template <class T>
inline Array<T, 3>& Array<T, 3>::operator~()
{
    apply(Complement<T>());
    return *this;
}

template <class T>
Array<T, 3>& Array<T, 3>::operator+=(const Array<T, 3>& a)
{
    Expects(extents == a.extents);
    for (size_type i = 0; i < size(); ++i) {
        elems[i] += a.data()[i];
    }
    return *this;
}

template <class T>
Array<T, 3>& Array<T, 3>::operator-=(const Array<T, 3>& a)
{
    Expects(extents == a.extents);
    for (size_type i = 0; i < size(); ++i) {
        elems[i] -= a.data()[i];
    }
    return *this;
}

template <class T>
inline Int_t Array<T, 3>::index(size_type i, size_type j, size_type k) const
{
    return i + j * strides[0] + k * strides[1];
}

template <class T>
void Array<T, 3>::assign(initializer_list_3d ilist)
{
    size_type n3 = ilist.size();
    size_type n1 = ilist.begin()->size();
    size_type n2 = ilist.begin()->begin()->size();

    elems.resize(n1 * n2 * n3);
    extents = {n1, n2, n3};
    strides = {n1, n1 * n2};

    size_type i = 0;
    size_type j = 0;
    size_type k = 0;

    for (const auto& iil : ilist) {
        for (const auto& il : iil) {
            for (const auto& v : il) {
                elems[index(i, j, k)] = v;
                ++j;
            }
            j = 0;
            ++i;
        }
        i = 0;
        ++k;
    }
}

}  // namespace srs

#endif  // SRS_ARRAY3_H
