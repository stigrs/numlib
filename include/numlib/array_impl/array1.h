// Copyright (c) 2018 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#ifndef NUMLIB_ARRAY1_H
#define NUMLIB_ARRAY1_H

#include <numlib/array_impl/array_ref.h>
#include <numlib/array_impl/functors.h>
#include <vector>
#include <cassert>

namespace Numlib {

// One-dimensional dense array (vector) class.
//
template <class T>
class Array<T, 1> {
public:
    static constexpr std::size_t rank = 1;

    using value_type = T;
    using size_type = Index;

    // Constructors:

    Array() : elems() {}

    explicit Array(size_type n) : elems(n) {}

    Array(size_type n, const T& value) : elems(n, value) {}

    Array(size_type n, T* ptr);

    // Copy elements referenced by Array_ref.
    template <class U>
    Array(const Array_ref<U, 1>& a);

    // Assignments:

    template <class U>
    Array& operator=(const Array_ref<U, 1>& a);

    // Element access:

    T& operator()(size_type i);
    const T& operator()(size_type i) const;

    T& operator[](size_type i);
    const T& operator[](size_type i) const;

    // Capacity:

    bool empty() const { return elems.empty(); }
    size_type size() const { return elems.size(); }

    // Modifiers:

    void clear() { elems.clear(); }

    void swap(Array& a) { elems.swap(a.elems); }

    void resize(size_type n) { elems.resize(n); }
    void resize(size_type n, const T& value) { elems.resize(n, value); }

    void push_back(const T& value) { elems.push_back(value); }
    void push_back(T&& value) { elems.push_back(value); }

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

    Array& operator+=(const Array& a);
    Array& operator-=(const Array& a);

private:
    std::vector<T> elems;  // storage
};

template <class T>
Array<T, 1>::Array(size_type n, T* ptr) : elems(n)
{
    for (size_type i = 0; i < size(); ++i) {
        elems[i] = ptr[i];
    }
}

template <class T>
template <class U>
Array<T, 1>::Array(const Array_ref<U, 1>& a) : elems(a.size())
{
    for (size_type i = 0; i < size(); ++i) {
        elems[i] = a[i];
    }
}

template <class T>
template <class U>
Array<T, 1>& Array<T, 1>::operator=(const Array_ref<U, 1>& a)
{
    resize(a.size());

    for (size_type i = 0; i < size(); ++i) {
        elems[i] = a[i];
    }
    return *this;
}

template <class T>
inline T& Array<T, 1>::operator()(size_type i)
{
#ifdef NDEBUG
    return elems[i];
#else
    return elems.at(i);
#endif
}

template <class T>
inline const T& Array<T, 1>::operator()(size_type i) const
{
#ifdef NDEBUG
    return elems[i];
#else
    return elems.at(i);
#endif
}

template <class T>
inline T& Array<T, 1>::operator[](size_type i)
{
#ifdef NDEBUG
    return elems[i];
#else
    return elems.at(i);
#endif
}

template <class T>
inline const T& Array<T, 1>::operator[](size_type i) const
{
#ifdef NDEBUG
    return elems[i];
#else
    return elems.at(i);
#endif
}

template <class T>
template <class F>
inline Array<T, 1>& Array<T, 1>::apply(F f)
{
    for (auto& v : elems) {
        f(v);
    }
    return *this;
}

template <class T>
template <class F>
inline Array<T, 1>& Array<T, 1>::apply(F f, const T& value)
{
    for (auto& v : elems) {
        f(v, value);
    }
    return *this;
}

template <class T>
inline Array<T, 1>& Array<T, 1>::operator=(const T& value)
{
    apply(Assign<T>(), value);
    return *this;
}

template <class T>
inline Array<T, 1>& Array<T, 1>::operator*=(const T& value)
{
    apply(Mul_assign<T>(), value);
    return *this;
}

template <class T>
inline Array<T, 1>& Array<T, 1>::operator/=(const T& value)
{
    apply(Div_assign<T>(), value);
    return *this;
}

template <class T>
inline Array<T, 1>& Array<T, 1>::operator%=(const T& value)
{
    apply(Mod_assign<T>(), value);
    return *this;
}

template <class T>
inline Array<T, 1>& Array<T, 1>::operator+=(const T& value)
{
    apply(Add_assign<T>(), value);
    return *this;
}

template <class T>
inline Array<T, 1>& Array<T, 1>::operator-=(const T& value)
{
    apply(Minus_assign<T>(), value);
    return *this;
}

template <class T>
inline Array<T, 1>& Array<T, 1>::operator&=(const T& value)
{
    apply(And_assign<T>(), value);
    return *this;
}

template <class T>
inline Array<T, 1>& Array<T, 1>::operator|=(const T& value)
{
    apply(Or_assign<T>(), value);
    return *this;
}

template <class T>
inline Array<T, 1>& Array<T, 1>::operator^=(const T& value)
{
    apply(Xor_assign<T>(), value);
    return *this;
}

template <class T>
Array<T, 1>& Array<T, 1>::operator+=(const Array<T, 1>& a)
{
    assert(size() == a.size());
    for (size_type i = 0; i < size(); ++i) {
        elems[i] += a[i];
    }
    return *this;
}

template <class T>
Array<T, 1>& Array<T, 1>::operator-=(const Array<T, 1>& a)
{
    assert(size() == a.size());
    for (size_type i = 0; i < size(); ++i) {
        elems[i] -= a[i];
    }
    return *this;
}

} // namespace Numlib

#endif // NUMLIB_ARRAY1_H

