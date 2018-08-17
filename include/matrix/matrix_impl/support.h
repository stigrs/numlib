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
#include <functional>
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

// Compute strides needed for subscript calculation and the number of elements
// given the extents.
//
// Storage order: Row-major
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

//
#endif  // NUMLIB_MATRIX_SUPPORT_H
