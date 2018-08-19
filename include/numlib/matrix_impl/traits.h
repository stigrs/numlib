// Copyright (c) 2008-2010 Kent State University
// Copyright (c) 2011-2012 Texas A&M University
// Copyright (c) Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#ifndef NUMLIB_MATRIX_TRAITS_H
#define NUMLIB_MATRIX_TRAITS_H

namespace Numlib {

// Deduce if the expression m.order is valid.
template <typename T>
struct get_order_result {
private:
    template <typename M>
    static auto check(const M& m) -> decltype(m.order);
    static substitution_failure check(...);

public:
    using type = decltype(check(std::declval<T>()));
};

template <typename T>
struct has_order : substitution_succeeded<typename get_order_result<T>::type> {
};

// Returns true if type is Matrix_type.
template <typename T>
constexpr bool Matrix_type()
{
    return has_order<T>::value;
}

namespace Matrix_impl {

    // Describes the structure of a nested std::initializer_list with
    // Matrix_init<T, N - 1> as its member type.
    template <typename T, std::size_t N>
    struct Matrix_init {
        using type =
            std::initializer_list<typename Matrix_init<T, N - 1>::type>;
    };

    // The N == 1 is special; that is were we go to the (most deeply nested)
    // std::initializer_list<T>.
    template <typename T>
    struct Matrix_init<T, 1> {
        using type = std::initializer_list<T>;
    };

    // To avoid surprises, N == 0 is defined to be an error.
    template <typename T>
    struct Matrix_init<T, 0>;

} // namespace Matrix_impl

} // namespace Numlib

#endif // NUMLIB_MATRIX_TRAITS_H
