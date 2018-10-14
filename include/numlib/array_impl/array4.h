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

#ifndef SRS_ARRAY4_H
#define SRS_ARRAY4_H

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
// Four-dimensional dense array class.
//
template <class T>
class Array<T, 4> {
public:
    static constexpr int rank = 4;

    typedef T value_type;
    typedef Int_t size_type;
    typedef typename std::vector<T>::iterator iterator;
    typedef typename std::vector<T>::const_iterator const_iterator;
    // clang-format off
    typedef typename std::initializer_list<std::initializer_list<std::initializer_list<std::initializer_list<T>>>> initializer_list_4d;
    // clang-format on

    // Constructors:

    Array() : elems(), extents{0, 0, 0, 0}, strides{0, 0, 0} {}

    Array(size_type n1, size_type n2, size_type n3, size_type n4)
        : elems(n1 * n2 * n3 * n4),
          extents{n1, n2, n3, n4},
          strides{n1, n1 * n2, n1 * n2 * n3}
    {
    }

    Array(
        size_type n1, size_type n2, size_type n3, size_type n4, const T& value)
        : elems(n1 * n2 * n3 * n4, value),
          extents{n1, n2, n3, n4},
          strides{n1, n1 * n2, n1 * n2 * n3}
    {
    }

    Array(size_type n1, size_type n2, size_type n3, size_type n4, T* ptr);

    template <Int_t n1, Int_t n2, Int_t n3, Int_t n4>
    Array(const T (&a)[n1][n2][n3][n4]);

    Array(initializer_list_4d ilist) { assign(ilist); }

    // T f(const T&) would be a typical type for f.
    template <class F>
    Array(const Array& a, F f);

    // T f(const T&, const Arg&) would be a typical type for f.
    template <class F, class Arg>
    Array(const Array& a, F f, const Arg& value);

    // Assignments:

    Array& operator=(initializer_list_4d ilist);

    // Element access:

    T& at(size_type i, size_type j, size_type k, size_type l);
    const T& at(size_type i, size_type j, size_type k, size_type l) const;

    T& operator()(size_type i, size_type j, size_type k, size_type l);
    const T& operator()(size_type i,
                        size_type j,
                        size_type k,
                        size_type l) const;

    Array_ref<T, 3> operator[](size_type i) { return slice(i); }
    Array_ref<const T, 3> operator[](size_type i) const { return slice(i); }

    // Iterators:

    iterator begin() { return elems.begin(); }
    iterator end() { return elems.end(); }

    const_iterator begin() const { return elems.begin(); }
    const_iterator end() const { return elems.end(); }

    // Slicing:

    Array_ref<T, 3> slice(size_type i);
    Array_ref<const T, 3> slice(size_type i) const;

    // Capacity:

    bool empty() const { return elems.empty(); }

    size_type size() const { return elems.size(); }
    size_type max_size() const { return elems.max_size(); }
    size_type capacity() const { return elems.capacity(); }
    size_type dim1() const { return extents[0]; }
    size_type dim2() const { return extents[1]; }
    size_type dim3() const { return extents[2]; }
    size_type dim4() const { return extents[3]; }
    size_type extent(size_type dim) const
    {
        Expects(dim >= 0 && dim < rank);
        return extents[dim];
    }

    // Modifiers:

    void clear();

    void swap(Array& a);

    void resize(size_type n1, size_type n2, size_type n3, size_type n4);
    void resize(
        size_type n1, size_type n2, size_type n3, size_type n4, const T& value);

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
    std::array<size_type, 4> extents;
    std::array<size_type, 3> strides;

    // Compute index.
    Int_t index(size_type i, size_type j, size_type k, size_type l) const;

    // Helper function for assigning initializer list.
    void assign(initializer_list_4d ilist);
};

template <class T>
Array<T, 4>::Array(
    size_type n1, size_type n2, size_type n3, size_type n4, T* ptr)
    : elems(n1 * n2 * n3 * n4),
      extents{n1, n2, n3, n4},
      strides{n1, n1 * n2, n1 * n2 * n3}
{
    for (size_type i = 0; i < size(); ++i) {
        elems[i] = ptr[i];
    }
}

template <class T>
template <Int_t n1, Int_t n2, Int_t n3, Int_t n4>
Array<T, 4>::Array(const T (&a)[n1][n2][n3][n4])
    : elems(n1 * n2 * n3 * n4),
      extents{n1, n2, n3, n4},
      strides{n1, n1 * n2, n1 * n2 * n3}
{
    for (size_type i = 0; i < extents[0]; ++i) {
        for (size_type j = 0; j < extents[1]; ++j) {
            for (size_type k = 0; k < extents[2]; ++k) {
                for (size_type l = 0; l < extents[3]; ++l) {
                    (*this)(i, j, k, l) = a[i][j][k][l];
                }
            }
        }
    }
}

template <class T>
template <class F>
Array<T, 4>::Array(const Array& a, F f)
{
    Expects(size() == a.size());
    for (size_type i = 0; i < size(); ++i) {
        elems[i] = f(a.data()[i]);
    }
}

template <class T>
template <class F, class Arg>
Array<T, 4>::Array(const Array& a, F f, const Arg& value)
{
    Expects(size() == a.size());
    for (size_type i = 0; i < size(); ++i) {
        elems[i] = f(a.data()[i], value);
    }
}

template <class T>
inline Array<T, 4>& Array<T, 4>::operator=(initializer_list_4d ilist)
{
    assign(ilist);
    return *this;
}

template <class T>
inline T& Array<T, 4>::at(size_type i, size_type j, size_type k, size_type l)
{
    Expects(i >= 0 && i < extents[0]);
    Expects(j >= 0 && j < extents[1]);
    Expects(k >= 0 && k < extents[2]);
    Expects(l >= 0 && l < extents[3]);
    return elems[index(i, j, k, l)];
}

template <class T>
inline const T& Array<T, 4>::at(size_type i,
                                size_type j,
                                size_type k,
                                size_type l) const
{
    Expects(i >= 0 && i < extents[0]);
    Expects(j >= 0 && j < extents[1]);
    Expects(k >= 0 && k < extents[2]);
    return elems[index(i, j, k, l)];
}

template <class T>
inline T& Array<T, 4>::operator()(size_type i,
                                  size_type j,
                                  size_type k,
                                  size_type l)
{
#ifdef NDEBUG
    return elems[index(i, j, k, l)];
#else
    return at(i, j, k, l);
#endif
}

template <class T>
inline const T& Array<T, 4>::operator()(size_type i,
                                        size_type j,
                                        size_type k,
                                        size_type l) const
{
#ifdef NDEBUG
    return elems[index(i, j, k, l)];
#else
    return at(i, j, k, l);
#endif
}

template <class T>
inline Array_ref<T, 3> Array<T, 4>::slice(size_type i)
{
    Expects(i >= 0 && i < extents[3]);
    return Array_ref<T, 3>(extents[0],
                           extents[1],
                           extents[2],
                           strides[0],
                           strides[1],
                           data() + i * strides[2]);
}

template <class T>
inline Array_ref<const T, 3> Array<T, 4>::slice(size_type i) const
{
    Expects(i >= 0 && i < extents[3]);
    return Array_ref<const T, 3>(extents[0],
                                 extents[1],
                                 extents[2],
                                 strides[0],
                                 strides[1],
                                 data() + i * strides[2]);
}

template <class T>
inline void Array<T, 4>::clear()
{
    elems.clear();
    extents = {0, 0, 0, 0};
    strides = {0, 0, 0};
}

template <class T>
inline void Array<T, 4>::swap(Array<T, 4>& a)
{
    elems.swap(a.elems);
    std::swap(extents, a.extents);
    std::swap(strides, a.strides);
}

template <class T>
inline void Array<T, 4>::resize(size_type n1,
                                size_type n2,
                                size_type n3,
                                size_type n4)
{
    elems.resize(n1 * n2 * n3 * n4);
    extents = {n1, n2, n3, n4};
    strides = {n1, n1 * n2, n1 * n2 * n3};
}

template <class T>
inline void Array<T, 4>::resize(
    size_type n1, size_type n2, size_type n3, size_type n4, const T& value)
{
    elems.resize(n1 * n2 * n3 * n4, value);
    extents = {n1, n2, n3, n4};
    strides = {n1, n1 * n2, n1 * n2 * n3};
}

template <class T>
template <class F>
inline Array<T, 4>& Array<T, 4>::apply(F f)
{
    for (auto& v : elems) {
        f(v);
    }
    return *this;
}

template <class T>
template <class F>
inline Array<T, 4>& Array<T, 4>::apply(F f, const T& value)
{
    for (auto& v : elems) {
        f(v, value);
    }
    return *this;
}

template <class T>
inline Array<T, 4>& Array<T, 4>::operator=(const T& value)
{
    apply(Assign<T>(), value);
    return *this;
}

template <class T>
inline Array<T, 4>& Array<T, 4>::operator*=(const T& value)
{
    apply(Mul_assign<T>(), value);
    return *this;
}

template <class T>
inline Array<T, 4>& Array<T, 4>::operator/=(const T& value)
{
    apply(Div_assign<T>(), value);
    return *this;
}

template <class T>
inline Array<T, 4>& Array<T, 4>::operator%=(const T& value)
{
    apply(Mod_assign<T>(), value);
    return *this;
}

template <class T>
inline Array<T, 4>& Array<T, 4>::operator+=(const T& value)
{
    apply(Add_assign<T>(), value);
    return *this;
}

template <class T>
inline Array<T, 4>& Array<T, 4>::operator-=(const T& value)
{
    apply(Minus_assign<T>(), value);
    return *this;
}

template <class T>
inline Array<T, 4>& Array<T, 4>::operator&=(const T& value)
{
    apply(And_assign<T>(), value);
    return *this;
}

template <class T>
inline Array<T, 4>& Array<T, 4>::operator|=(const T& value)
{
    apply(Or_assign<T>(), value);
    return *this;
}

template <class T>
inline Array<T, 4>& Array<T, 4>::operator^=(const T& value)
{
    apply(Xor_assign<T>(), value);
    return *this;
}

template <class T>
inline Array<T, 4>& Array<T, 4>::operator!()
{
    apply(Not<T>());
    return *this;
}

template <class T>
inline Array<T, 4>& Array<T, 4>::operator-()
{
    apply(Unary_minus<T>());
    return *this;
}

template <class T>
inline Array<T, 4>& Array<T, 4>::operator~()
{
    apply(Complement<T>());
    return *this;
}

template <class T>
Array<T, 4>& Array<T, 4>::operator+=(const Array<T, 4>& a)
{
    Expects(extents == a.extents);
    for (size_type i = 0; i < size(); ++i) {
        elems[i] += a.data()[i];
    }
    return *this;
}

template <class T>
Array<T, 4>& Array<T, 4>::operator-=(const Array<T, 4>& a)
{
    Expects(extents == a.extents);
    for (size_type i = 0; i < size(); ++i) {
        elems[i] -= a.data()[i];
    }
    return *this;
}

template <class T>
inline Int_t Array<T, 4>::index(size_type i,
                                size_type j,
                                size_type k,
                                size_type l) const
{
    return i + j * strides[0] + k * strides[1] + l * strides[2];
}

template <class T>
void Array<T, 4>::assign(initializer_list_4d ilist)
{
    size_type n4 = ilist.size();
    size_type n3 = ilist.begin()->size();
    size_type n1 = ilist.begin()->begin()->size();
    size_type n2 = ilist.begin()->begin()->begin()->size();

    elems.resize(n1 * n2 * n3 * n4);
    extents = {n1, n2, n3, n4};
    strides = {n1, n1 * n2, n1 * n2 * n3};

    size_type i = 0;
    size_type j = 0;
    size_type k = 0;
    size_type l = 0;

    for (const auto& iiil : ilist) {
        for (const auto& iil : iiil) {
            for (const auto& il : iil) {
                for (const auto& v : il) {
                    elems[index(i, j, k, l)] = v;
                    ++j;
                }
                j = 0;
                ++i;
            }
            i = 0;
            ++k;
        }
        k = 0;
        ++l;
    }
}

}  // namespace srs

#endif  // SRS_ARRAY4_H
