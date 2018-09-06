// Copyright (c) 2008-2010 Kent State University
// Copyright (c) 2011-2012 Texas A&M University
// Copyright (c) 2018 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#ifndef NUMLIB_TRAITS_H
#define NUMLIB_TRAITS_H

#include <cstddef>
#include <type_traits>
#include <initializer_list>
#include <utility>

//------------------------------------------------------------------------------

// Signed index type:
//
// C++ Core Guidelines recommends to use a signed integer type for subscripts/
// indices. The Guidelines Support Library (GSL) provides gsl::index as a
// typedef to std::ptrdiff_t to avoid the ugliness of using ptrdiff_t.
//
// Here, a typedef to ptrdiff_t is provided in the global namespace in order to
// avoid linking with the GSL. Hopefully, this can be replaced by std::index in
// the future.

using Index = std::ptrdiff_t;

//------------------------------------------------------------------------------

// Type queries:

template <typename U>
using Value_type = typename U::value_type;

template <typename T>
constexpr bool Integer_type()
{
    return std::is_integral<T>::value;
}

template <typename T>
constexpr bool Real_type()
{
    return std::is_floating_point<T>::value;
}

//------------------------------------------------------------------------------

// An alias to U if T has type const U, otherwise T.
template <typename T>
using Remove_const = typename std::remove_const<T>::type;

//------------------------------------------------------------------------------

// Enable if B is true.
template <bool B, typename T = void>
using Enable_if = typename std::enable_if<B, T>::type;

//------------------------------------------------------------------------------

// Return true if T is convertible to U.
template <typename T, typename U>
constexpr bool Convertible()
{
    return std::is_convertible<T, U>::value;
}

//------------------------------------------------------------------------------

// Type predicates:

// Return true if every argument is true of if no arguments are given.

constexpr bool All() { return true; }

template <typename... Args>
constexpr bool All(bool b, Args... args)
{
    return b && All(args...);
}

// Return true if some (at least one) argument is true.

constexpr bool Some() { return false; }

template <typename... Args>
constexpr bool Some(bool b, Args... args)
{
    return b || Some(args...);
}

// Return true when types T and U are the same type.
template <typename T, typename U>
constexpr bool Same()
{
    return std::is_same<T, U>::value;
}

//------------------------------------------------------------------------------

// SFINAE support:

struct substitution_failure { // represent a failure to declare something
};

template <typename T>
struct substitution_succeeded : std::true_type {
};

template <>
struct substitution_succeeded<substitution_failure> : std::false_type {
};

//------------------------------------------------------------------------------
//
// Type cast:

// A searchable way to do narrowing casts of values.
template <typename T, typename U>
constexpr T narrow_cast(U&& u)
{
    return static_cast<T>(std::forward<U>(u));
}

#endif // NUMLIB_TRAITS_H
