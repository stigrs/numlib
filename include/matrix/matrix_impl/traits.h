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

#ifndef NUMLIB_MATRIX_TRAITS_H
#define NUMLIB_MATRIX_TRAITS_H

#include <type_traits>

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

// Return true if every argument is true of if no arguments are given.

constexpr bool All() { return true; }

template <typename... Args>
constexpr bool All(bool b, Args... args)
{
    return b && All(args...);
}

#endif  // NUMLIB_MATRIX_TRAITS_H
