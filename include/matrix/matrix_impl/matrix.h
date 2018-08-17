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

#ifndef NUMLIB_MATRIX_MATRIX_H
#define NUMLIB_MATRIX_MATRIX_H

#include <matrix/matrix_impl/matrix_base.h>
#include <matrix/matrix_impl/support.h>
#include <initializer_list>
#include <iostream>
#include <vector>


namespace numlib {

//------------------------------------------------------------------------------

template <typename T, std::size_t N>
using Matrix_initializer = typename matrix_impl::Matrix_init<T, N>::type;

//------------------------------------------------------------------------------

template <typename T, std::size_t N>
class Matrix : public Matrix_base<T, N> {
public:
    using size_type      = typename Matrix_base<T, N>::size_type;
    using iterator       = typename std::vector<T>::iterator;
    using const_iterator = typename std::vector<T>::const_iterator;

    Matrix()  = default;
    ~Matrix() = default;

    // Move construction and assignment:
    Matrix(Matrix&&) = default;
    Matrix& operator=(Matrix&&) = default;

    // Copy construction and assignment:
    Matrix(const Matrix&) = default;
    Matrix& operator=(const Matrix&) = default;

    // Specify the extents:
    template <typename... Exts>
    explicit Matrix(Exts... exts);

    // Initialize and assign from list:

    template <typename U>
    Matrix(std::initializer_list<U>) = delete;  // don't use {} except for elems

    template <typename U>
    Matrix& operator=(std::initializer_list<U>) = delete;

    Matrix(Matrix_initializer<T, N>);
    Matrix& operator=(Matrix_initializer<T, N>);

    // Total number of elements.
    size_type size() const { return elems.size(); }

    // "Flat" element access:
    T* data() { return elems.data(); }
    const T* data() const { return elems.data(); }

    // clang-format off
    template<typename... Args>
    Enable_if<matrix_impl::Requesting_element<Args...>(), T&>
    operator()(Args... args)
    {
        return Matrix_base<T, N>::template operator()<Args...>(args...);
    }

    template<typename... Args>
    Enable_if<matrix_impl::Requesting_element<Args...>(), const T&>
    operator()(Args... args) const
    {
        return Matrix_base<T, N>::template operator()<Args...>(args...);
    }
    // clang-format on
private:
    std::vector<T> elems;
};

template <typename T, std::size_t N>
template <typename... Exts>
Matrix<T, N>::Matrix(Exts... exts)
    : Matrix_base<T, N>{exts...}, elems(this->desc.size)
{
}

template <typename T, std::size_t N>
Matrix<T, N>::Matrix(Matrix_initializer<T, N> init)
{
    this->desc.extents = matrix_impl::derive_extents<N>(init);
    matrix_impl::compute_strides(this->desc);
    elems.reserve(this->desc.size);
    matrix_impl::insert_flat(init, elems);
    assert(elems.size() == this->desc.size);
}

}  // namespace numlib

#endif  // NUMLIB_MATRIX_MATRIX_H
