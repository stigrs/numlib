////////////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2018 Stig Rune Sellevag. All rights reserved.
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

#ifndef NUMLIB_TRAITS_H
#define NUMLIB_TRAITS_H

#include <numlib/matrix_impl/matrix_fwd.h>
#include <type_traits>

//------------------------------------------------------------------------------

template <typename U>
using Value_type = typename U::value_type;

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

struct substitution_failure {  // represent a failure to declare something
};

template <typename T>
struct substitution_succeeded : std::true_type {
};

template <>
struct substitution_succeeded<substitution_failure> : std::false_type {
};

// Return true if type is Matrix_type:

template <typename M>
struct get_matrix_type_result {
    template <typename T, std::size_t N, typename = Enable_if<(N >= 1)>>
    static bool check(const Numlib::Matrix<T, N>& m);

    template <typename T, std::size_t N, typename = Enable_if<(N >= 1)>>
    static bool check(const Numlib::Matrix_ref<T, N>& m);

    static substitution_failure check(...);

    using type = decltype(check(std::declval<M>()));
};

template <typename T>
struct has_matrix_type
    : substitution_succeeded<typename get_matrix_type_result<T>::type> {
};

template <typename M>
constexpr bool Matrix_type()
{
    return has_matrix_type<M>::value;
}

#endif  // NUMLIB_TRAITS_H
