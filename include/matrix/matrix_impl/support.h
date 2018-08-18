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

#ifndef NUMLIB_MATRIX_SUPPORT_H
#define NUMLIB_MATRIX_SUPPORT_H

#include <matrix/matrix_impl/traits.h>
#include <algorithm>
#include <array>
#include <cassert>
#include <functional>
#include <initializer_list>
#include <numeric>

//------------------------------------------------------------------------------

// Forward declarations:

template <std::size_t N>
struct Matrix_slice;

//------------------------------------------------------------------------------

namespace matrix_impl {

template <typename... Args>
constexpr bool Requesting_element()
{
    return All(Convertible<Args, std::size_t>()...);
}

//------------------------------------------------------------------------------

// Matrix list initialization:

// Forward declaration:

template <std::size_t N, typename List>
bool check_non_jagged(const List& list);

// Describes the structure of a nested std::initializer_list with
// Matrix_init<T, N - 1> as its member type.
template <typename T, std::size_t N>
struct Matrix_init {
    using type = std::initializer_list<typename Matrix_init<T, N - 1>::type>;
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

template <std::size_t N, typename I, typename List>
Enable_if<(N == 1), void> add_extents(I& first, const List& list)
{
    *first++ = list.size();
}

// Recursion through nested std::initializer_list.
template <std::size_t N, typename I, typename List>
Enable_if<(N > 1), void> add_extents(I& first, const List& list)
{
    assert(check_non_jagged<N>(list));
    *first++ = list.size();  // store this size (extent)
    add_extents<N - 1>(first, *list.begin());
}

// Determine the shape of the Matrix:
//   + Checks that the tree is really N deep
//   + Checks that each row has the same number of elements
//   + Sets the extent of each row
//
template <std::size_t N, typename List>
std::array<std::size_t, N> derive_extents(const List& list)
{
    std::array<std::size_t, N> a;
    auto f = a.begin();
    add_extents<N>(f, list);  // add sizes (extents) to a
    return a;
}

// Check that all rows have the same number of elements.
template <std::size_t N, typename List>
bool check_non_jagged(const List& list)
{
    auto i = list.begin();
    for (auto j = i + 1; j != list.end(); ++j) {
        if (derive_extents<N - 1>(*i) != derive_extents<N - 1>(*j)) {
            return false;
        }
    }
    return true;
}


// When we reach a list with non-initializer_list elements, we insert
// those elements into our vector.
template <typename T, typename Vec>
void add_list(const T* first, const T* last, Vec& vec)
{
    vec.insert(vec.end(), first, last);
}

template <typename T, typename Vec>
void add_list(const std::initializer_list<T>* first,
              const std::initializer_list<T>* last,
              Vec& vec)
{
    while (first != last) {
        add_list(first->begin(), first->end(), vec);
        ++first;
    }
}

// Copy elements of the tree of std::initializer_list to a Matrix<T, N>.
template <typename T, typename Vec>
void insert_flat(std::initializer_list<T> list, Vec& vec)
{
    add_list(list.begin(), list.end(), vec);
}

//------------------------------------------------------------------------------

// Compute strides needed for subscript calculation and the number of elements
// given the extents.
//
// Note: Row-major storage order.
//
template <std::size_t N>
void compute_strides(Matrix_slice<N>& ms)
{
    ms.strides[N - 1] = 1;                      // last stride is 1
    for (std::size_t i = N - 1; i != 0; --i) {  // compute stride for each dim
        ms.strides[i - 1] = ms.strides[i] * ms.extents[i];
    }
    ms.size = ms.extents[0] * ms.strides[0];
}

// Compute total number of elements given the extents.
template <std::size_t N>
std::size_t compute_size(const std::array<std::size_t, N>& exts)
{
    return std::accumulate(
        exts.begin(), exts.end(), 1, std::multiplies<std::size_t>{});
}

// Return true if each element in range is within the bounds of the
// corresponding extent.
template <std::size_t N, typename... Dims>
bool check_bounds(const Matrix_slice<N>& slice, Dims... dims)
{
    std::size_t indexes[N]{std::size_t(dims)...};
    return std::equal(
        indexes, indexes + N, slice.extents.begin(), std::less<std::size_t>{});
}

}  // namespace matrix_impl

#endif  // NUMLIB_MATRIX_SUPPORT_H
