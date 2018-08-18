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

#ifndef NUMLIB_MATRIX_MATRIX_H
#define NUMLIB_MATRIX_MATRIX_H

#include <numlib/matrix_impl/matrix_base.h>
#include <numlib/matrix_impl/support.h>
#include <initializer_list>
#include <utility>
#include <vector>


namespace Numlib {

//------------------------------------------------------------------------------

template <typename T, std::size_t N>
using Matrix_initializer = typename Matrix_impl::Matrix_init<T, N>::type;

//------------------------------------------------------------------------------

template <typename T, std::size_t N>
class Matrix : public Matrix_base<T, N> {
public:
    using size_type      = typename Matrix_base<T, N>::size_type;
    using iterator       = typename std::vector<T>::iterator;
    using const_iterator = typename std::vector<T>::const_iterator;

    Matrix() = default;

    // Move construction and assignment:
    Matrix(Matrix&&) = default;
    Matrix& operator=(Matrix&&) = default;

    // Copy construction and assignment:
    Matrix(const Matrix&) = default;
    Matrix& operator=(const Matrix&) = default;

    // Specify the extents:
    template <typename... Exts>
    explicit Matrix(Exts... exts);

// Construct and assign from Matrix_ref:
#if 0
    template <typename U>
    Matrix(const Matrix_ref<U, N>&);

    template <typename U>
    Matrix& operator=(const Matrix_ref<U, N>&);
#endif
    // Initialize and assign from list:

    template <typename U>
    Matrix(std::initializer_list<U>) = delete;  // don't use {} except for elems

    template <typename U>
    Matrix& operator=(std::initializer_list<U>) = delete;

    Matrix(Matrix_initializer<T, N>);
    Matrix& operator=(Matrix_initializer<T, N>);

    ~Matrix() = default;

    // "Flat" element access:
    T* data() { return elems.data(); }
    const T* data() const { return elems.data(); }

    // Properties:

    bool empty() const { return elems.empty(); }

    // Subscripting:

    // clang-format off
    template <typename... Args>
	Enable_if<Matrix_impl::Requesting_element<Args...>(), T&> 
	operator()(Args... args)
	{
		assert(Matrix_impl::check_bounds(this->desc, args...));
		return *(data() + this->desc(args...));
	}

    template <typename... Args>
	Enable_if<Matrix_impl::Requesting_element<Args...>(), const T&> 
	operator()(Args... args) const 
	{
		assert(Matrix_impl::check_bounds(this->desc, args...));
		return *(data() + this->desc(args...));
	}
    // clang-format on

    // Iterators:

    iterator begin() { return elems.begin(); }
    iterator end() { return elems.end(); }

    const_iterator begin() const { return elems.begin(); }
    const_iterator end() const { return elems.end(); }

    // Mutators:

    void swap(Matrix& m);

    // Apply f(x) for every element x:
    template <typename F>
    Matrix& apply(F f);

    // Apply f(x, mx) for corresponding elements *this and m:
    template <typename M, typename F>
    Matrix& apply(const M& m, F f);

private:
    std::vector<T> elems;
};

template <typename T, std::size_t N>
template <typename... Exts>
inline Matrix<T, N>::Matrix(Exts... exts)
    : Matrix_base<T, N>{exts...}, elems(this->desc.size)
{
}

template <typename T, std::size_t N>
inline Matrix<T, N>::Matrix(Matrix_initializer<T, N> init)
{
    this->desc.start   = 0;
    this->desc.extents = Matrix_impl::derive_extents<N>(init);
    Matrix_impl::compute_strides(this->desc);
    elems.reserve(this->desc.size);
    Matrix_impl::insert_flat(init, elems);
    assert(elems.size() == this->desc.size);
}

template <typename T, std::size_t N>
inline Matrix<T, N>& Matrix<T, N>::operator=(Matrix_initializer<T, N> init)
{
    Matrix tmp(init);
    swap(tmp);
    return *this;
}

template <typename T, std::size_t N>
inline void Matrix<T, N>::swap(Matrix& m)
{
    std::swap(this->desc, m.desc);
    elems.swap(m.elems);
}

template <typename T, std::size_t N>
template <typename F>
Matrix<T, N>& Matrix<T, N>::apply(F f)
{
    for (auto& x : elems) {
        f(x);
    }
    return *this;
}

template <typename T, std::size_t N>
template <typename M, typename F>
Matrix<T, N>& Matrix<T, N>::apply(const M& m, F f)
{
    assert(same_extents(this->desc, m.descriptor()));
    for (auto i = begin(), j = m.begin(); i != end(); ++i, ++j) {
        f(*i, *j);
    }
    return *this;
}

}  // namespace Numlib

#endif  // NUMLIB_MATRIX_MATRIX_H
