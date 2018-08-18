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

#ifndef NUMLIB_MATRIX_MATRIX_BASE_H
#define NUMLIB_MATRIX_MATRIX_BASE_H

#include <numlib/matrix_impl/matrix_slice.h>
#include <numlib/matrix_impl/support.h>
#include <numlib/traits/traits.h>
#include <cassert>

namespace Numlib {

template <typename T, std::size_t N>
class Matrix_base {
public:
    static constexpr std::size_t rank = N;

    using value_type = T;
    using size_type  = std::size_t;

    Matrix_base() = default;

    template <typename... Exts>
    explicit Matrix_base(Exts... exts) : desc{exts...}
    {
    }

    explicit Matrix_base(const Matrix_slice<N>& s) : desc{s} {}

    // Move construction and assignment:
    Matrix_base(Matrix_base&&) = default;
    Matrix_base& operator=(Matrix_base&&) = default;

    // Copy construction and assignment:
    Matrix_base(const Matrix_base&) = default;
    Matrix_base& operator=(const Matrix_base&) = default;

    ~Matrix_base() = default;

    // Number of dimensions:
    static constexpr size_type order() { return rank; }

    // Total number of elements:
    size_type size() const { return desc.size; }

    // Number of elements in the N'th dimension:
    size_type extent(size_type n) const
    {
        assert(n >= 0 && n < rank);
        return desc.extents[n];
    }

    // Number of rows.
    size_type rows() const
    {
        assert(N >= 1);
        return desc.extents[0];
    }

    // Number of columns.
    size_type cols() const
    {
        assert(N >= 2);
        return desc.extents[1];
    }

    // The slice defining subscripting:
    const Matrix_slice<N>& descriptor() const { return desc; }

protected:
    Matrix_slice<N> desc;
};

}  // namespace Numlib

#endif  // NUMLIB_MATRIX_MATRIX_BASE_H
