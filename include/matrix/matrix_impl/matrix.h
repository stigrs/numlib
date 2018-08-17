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
#include <vector>

namespace numlib {

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

    // Total number of elements.
    size_type size() const { return elems.size(); }

    // "Flat" element access:
    T* data() { return elems.data(); }
    const T* data() const { return elems.data(); }

private:
    std::vector<T> elems;
};

template <typename T, std::size_t N>
template <typename... Exts>
Matrix<T, N>::Matrix(Exts... exts)
    : Matrix_base<T, N>{exts...}, elems(this->desc.size)
{
}

}  // namespace numlib

#endif  // NUMLIB_MATRIX_MATRIX_H
