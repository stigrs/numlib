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

#ifndef SRS_ARRAY_OPR_H
#define SRS_ARRAY_OPR_H

#include <algorithm>
#include <gsl/gsl>


namespace srs {

// Comparison operators:

template <class T, int N>
inline bool operator==(const Array<T, N>& a, const Array<T, N>& b)
{
    return std::equal(a.begin(), a.end(), b.begin());
}

template <class T, int N>
inline bool operator!=(const Array<T, N>& a, const Array<T, N>& b)
{
    return !(a == b);
}

template <class T, int N>
inline bool operator<(const Array<T, N>& a, const Array<T, N>& b)
{
    return std::lexicographical_compare(a.begin(), a.end(), b.begin(), b.end());
}

template <class T, int N>
inline bool operator>(const Array<T, N>& a, const Array<T, N>& b)
{
    return (b < a);
}

template <class T, int N>
inline bool operator<=(const Array<T, N>& a, const Array<T, N>& b)
{
    return !(a > b);
}

template <class T, int N>
inline bool operator>=(const Array<T, N>& a, const Array<T, N>& b)
{
    return !(a < b);
}

//------------------------------------------------------------------------------

// Array addition:

template <class T, int N>
inline Array<T, N> operator+(const Array<T, N>& a, const Array<T, N>& b)
{
    Array<T, N> result(a);
    return result += b;
}

template <class T, int N>
inline Array<T, N> operator+(const Array<T, N>& a, const Array_ref<T, N>& b)
{
    Array<T, N> result(a);
    return result += b;
}

template <class T, int N>
inline Array<T, N> operator+(const Array_ref<T, N>& a, const Array<T, N>& b)
{
    Array<T, N> result(a);
    return result += b;
}

template <class T, int N>
inline Array<T, N> operator+(const Array_ref<T, N>& a, const Array_ref<T, N>& b)
{
    Array<T, N> result(a);
    return result += b;
}

template <class T, int N>
inline Array<T, N> operator+(const Array_ref<const T, N>& a,
                             const Array_ref<const T, N>& b)
{
    Array<T, N> result(a);
    return result += b;
}

//------------------------------------------------------------------------------

// Array subtraction:

template <class T, int N>
inline Array<T, N> operator-(const Array<T, N>& a, const Array<T, N>& b)
{
    Array<T, N> result(a);
    return result -= b;
}

template <class T, int N>
inline Array<T, N> operator-(const Array<T, N>& a, const Array_ref<T, N>& b)
{
    Array<T, N> result(a);
    return result -= b;
}

template <class T, int N>
inline Array<T, N> operator-(const Array_ref<T, N>& a, const Array<T, N>& b)
{
    Array<T, N> result(a);
    return result -= b;
}

template <class T, int N>
inline Array<T, N> operator-(const Array_ref<T, N>& a, const Array_ref<T, N>& b)
{
    Array<T, N> result(a);
    return result -= b;
}

template <class T, int N>
inline Array<T, N> operator-(const Array_ref<const T, N>& a,
                             const Array_ref<const T, N>& b)
{
    Array<T, N> result(a);
    return result -= b;
}

//------------------------------------------------------------------------------

// Scalar addition:

template <class T, int N>
inline Array<T, N> operator+(const Array<T, N>& a, const T& scalar)
{
    Array<T, N> result(a);
    return result += scalar;
}

template <class T, int N>
inline Array<T, N> operator+(const T& scalar, const Array<T, N>& a)
{
    Array<T, N> result(a);
    return result += scalar;
}

template <class T, int N>
inline Array<T, N> operator+(const Array_ref<T, N>& a, const T& scalar)
{
    Array<T, N> result(a);
    return result += scalar;
}

template <class T, int N>
inline Array<T, N> operator+(const T& scalar, const Array_ref<T, N>& a)
{
    Array<T, N> result(a);
    return result += scalar;
}

template <class T, int N>
inline Array<T, N> operator+(const Array_ref<const T, N>& a, const T& scalar)
{
    Array<T, N> result(a);
    return result += scalar;
}

template <class T, int N>
inline Array<T, N> operator+(const T& scalar, const Array_ref<const T, N>& a)
{
    Array<T, N> result(a);
    return result += scalar;
}

//------------------------------------------------------------------------------

// Scalar subtraction:

template <class T, int N>
inline Array<T, N> operator-(const Array<T, N>& a, const T& scalar)
{
    Array<T, N> result(a);
    return result -= scalar;
}

template <class T, int N>
inline Array<T, N> operator-(const T& scalar, const Array<T, N>& a)
{
    Array<T, N> result(a);
    return result -= scalar;
}

template <class T, int N>
inline Array<T, N> operator-(const Array_ref<T, N>& a, const T& scalar)
{
    Array<T, N> result(a);
    return result -= scalar;
}

template <class T, int N>
inline Array<T, N> operator-(const T& scalar, const Array_ref<T, N>& a)
{
    Array<T, N> result(a);
    return result -= scalar;
}

template <class T, int N>
inline Array<T, N> operator-(const Array_ref<const T, N>& a, const T& scalar)
{
    Array<T, N> result(a);
    return result -= scalar;
}

template <class T, int N>
inline Array<T, N> operator-(const T& scalar, const Array_ref<const T, N>& a)
{
    Array<T, N> result(a);
    return result -= scalar;
}

//------------------------------------------------------------------------------

// Scalar multiplication:

template <class T, int N>
inline Array<T, N> operator*(const Array<T, N>& a, const T& scalar)
{
    Array<T, N> result(a);
    return result *= scalar;
}

template <class T, int N>
inline Array<T, N> operator*(const T& scalar, const Array<T, N>& a)
{
    Array<T, N> result(a);
    return result *= scalar;
}

template <class T, int N>
inline Array<T, N> operator*(const Array_ref<T, N>& a, const T& scalar)
{
    Array<T, N> result(a);
    return result *= scalar;
}

template <class T, int N>
inline Array<T, N> operator*(const T& scalar, const Array_ref<T, N>& a)
{
    Array<T, N> result(a);
    return result *= scalar;
}

template <class T, int N>
inline Array<T, N> operator*(const Array_ref<const T, N>& a, const T& scalar)
{
    Array<T, N> result(a);
    return result *= scalar;
}

template <class T, int N>
inline Array<T, N> operator*(const T& scalar, const Array_ref<const T, N>& a)
{
    Array<T, N> result(a);
    return result *= scalar;
}

//------------------------------------------------------------------------------

// Matrix-matrix multiplication:

// Declaration.
template <class A1, class A2, class A3>
void mm_mul(const A1& a, const A2& b, A3& c);

template <class T>
inline Array<T, 2> operator*(const Array<T, 2>& a, const Array<T, 2>& b)
{
    Array<T, 2> result;
    mm_mul(a, b, result);
    return result;
}

template <class T>
inline Array<T, 2> operator*(const Array_ref<T, 2>& a, const Array_ref<T, 2>& b)
{
    Array<T, 2> result;
    mm_mul(a, b, result);
    return result;
}

template <class T>
inline Array<T, 2> operator*(const Array_ref<const T, 2>& a,
                             const Array_ref<const T, 2>& b)
{
    Array<T, 2> result;
    mm_mul(a, b, result);
    return result;
}

template <class T>
inline Array<T, 2> operator*(const Array<T, 2>& a, const Array_ref<T, 2>& b)
{
    Array<T, 2> result;
    mm_mul(a, b, result);
    return result;
}

template <class T>
inline Array<T, 2> operator*(const Array<T, 2>& a,
                             const Array_ref<const T, 2>& b)
{
    Array<T, 2> result;
    mm_mul(a, b, result);
    return result;
}

template <class T>
inline Array<T, 2> operator*(const Array_ref<T, 2>& a, const Array<T, 2>& b)
{
    Array<T, 2> result;
    mm_mul(a, b, result);
    return result;
}

template <class T>
inline Array<T, 2> operator*(const Array_ref<const T, 2>& a,
                             const Array<T, 2>& b)
{
    Array<T, 2> result;
    mm_mul(a, b, result);
    return result;
}

//------------------------------------------------------------------------------

// Matrix-vector multiplication:

// Declaration.
template <class A1, class A2, class A3>
void mv_mul(const A1& a, const A2& v, A3& w);

template <class T>
inline Array<T, 1> operator*(const Array<T, 2>& a, const Array<T, 1>& v)
{
    Array<T, 1> result;
    mv_mul(a, v, result);
    return result;
}

template <class T>
inline Array<T, 1> operator*(const Array_ref<T, 2>& a, const Array_ref<T, 1>& v)
{
    Array<T, 1> result;
    mv_mul(a, v, result);
    return result;
}

template <class T>
inline Array<T, 1> operator*(const Array_ref<const T, 2>& a,
                             const Array_ref<const T, 1>& v)
{
    Array<T, 1> result;
    mv_mul(a, v, result);
    return result;
}

template <class T>
inline Array<T, 1> operator*(const Array<T, 2>& a, const Array_ref<T, 1>& v)
{
    Array<T, 1> result;
    mv_mul(a, v, result);
    return result;
}

template <class T>
inline Array<T, 1> operator*(const Array<T, 2>& a,
                             const Array_ref<const T, 1>& v)
{
    Array<T, 1> result;
    mv_mul(a, v, result);
    return result;
}

template <class T>
inline Array<T, 1> operator*(const Array_ref<T, 2>& a, const Array<T, 1>& v)
{
    Array<T, 1> result;
    mv_mul(a, v, result);
    return result;
}

template <class T>
inline Array<T, 1> operator*(const Array_ref<const T, 2>& a,
                             const Array<T, 1>& v)
{
    Array<T, 1> result;
    mv_mul(a, v, result);
    return result;
}

//------------------------------------------------------------------------------

// Matrix-matrix multiplication.
template <class A1, class A2, class A3>
// requires A1 = Array<T, 2>, A2 = Array<T, 2>, A3 = Array<T, 2>
void mm_mul(const A1& a, const A2& b, A3& c)
{
    using value_type = typename A1::value_type;
    using size_type  = typename A1::size_type;

    Expects(A1::rank == 2);
    Expects(A2::rank == 2);
    Expects(A3::rank == 2);
    Expects(a.cols() == b.rows());

    c.resize(a.rows(), b.cols());

    for (size_type j = 0; j < b.cols(); ++j) {
        for (size_type i = 0; i < a.rows(); ++i) {
            c(i, j) = value_type(0);
            for (size_type k = 0; k < a.cols(); ++k) {
                c(i, j) += a(i, k) * b(k, j);
            }
        }
    }
}

// Matrix-vector multiplication.
template <class A1, class A2, class A3>
// requires A1 = Array<T, 2>, A2 = Array<T, 1>, A3 = Array<T, 1>
void mv_mul(const A1& a, const A2& v, A3& w)
{
    using value_type = typename A1::value_type;
    using size_type  = typename A1::size_type;

    Expects(A1::rank == 2);
    Expects(A2::rank == 1);
    Expects(A3::rank == 1);
    Expects(v.size() == a.cols());

    w.resize(a.rows());
    w = value_type(0);

    for (size_type j = 0; j < a.cols(); ++j) {
        for (size_type i = 0; i < a.rows(); ++i) {
            w(i) += a(i, j) * v(j);
        }
    }
}

//------------------------------------------------------------------------------

// Algorithms:

// Swap arrays.
template <class T, int N>
inline void swap(Array<T, N>& a, Array<T, N>& b)
{
    a.swap(b);
}

// Sort vector.
template <class T>
inline void sort(Array<T, 1>& vec, bool ascending = true)
{
    if (ascending) {
        std::sort(vec.begin(), vec.end(), std::less<T>());
    }
    else {  // descending
        std::sort(vec.begin(), vec.end(), std::greater<T>());
    }
}

// Sort matrix.
template <typename T>
void sort(Array<T, 2>& a, int dim = 2, bool ascending = true)
{
    using size_type = typename Array<T, 2>::size_type;

    if (dim == 1) {  // sort elements along each row
        if (ascending) {
            for (size_type i = 0; i < a.rows(); ++i) {
                // This is not terribly elegant, but I don't know (yet) how
                // to make a random access slice iterator.
                Array<T, 1> ri = a.row(i);
                std::sort(ri.begin(), ri.end(), std::less<T>());
                a.row(i) = ri;
            }
        }
        else {
            for (size_type i = 0; i < a.rows(); ++i) {
                Array<T, 1> ri = a.row(i);
                std::sort(ri.begin(), ri.end(), std::greater<T>());
                a.row(i) = ri;
            }
        }
    }
    else {  // sort elements along each column
        if (ascending) {
            for (size_type j = 0; j < a.cols(); ++j) {
                Array<T, 1> cj = a.column(j);
                std::sort(cj.begin(), cj.end(), std::less<T>());
                a.column(j) = cj;
            }
        }
        else {
            for (size_type j = 0; j < a.cols(); ++j) {
                Array<T, 1> cj = a.column(j);
                std::sort(cj.begin(), cj.end(), std::greater<T>());
                a.column(j) = cj;
            }
        }
    }
}

// Copy n elements from a to b.
template <class T>
void copy_n(srs::size_t n, const Array<T, 1>& a, Array<T, 1>& b)
{
    if (n > 0) {
        for (srs::size_t i = 0; i < n; ++i) {
            b(i) = a(i);
        }
    }
}

template <class T>
void copy_n(srs::size_t n, const Array_ref<T, 1>& a, Array_ref<T, 1>& b)
{
    if (n > 0) {
        for (srs::size_t i = 0; i < n; ++i) {
            b(i) = a(i);
        }
    }
}

template <class T>
void copy_n(srs::size_t n, const Array_ref<const T, 1>& a, Array_ref<T, 1>& b)
{
    if (n > 0) {
        for (srs::size_t i = 0; i < n; ++i) {
            b(i) = a(i);
        }
    }
}

template <class T>
void copy_n(srs::size_t n, const Array<T, 1>& a, Array_ref<T, 1>& b)
{
    if (n > 0) {
        for (srs::size_t i = 0; i < n; ++i) {
            b(i) = a(i);
        }
    }
}

template <class T>
void copy_n(srs::size_t n, const Array_ref<T, 1>& a, Array<T, 1>& b)
{
    if (n > 0) {
        for (srs::size_t i = 0; i < n; ++i) {
            b(i) = a(i);
        }
    }
}

template <class T>
void copy_n(srs::size_t n, const Array_ref<const T, 1>& a, Array<T, 1>& b)
{
    if (n > 0) {
        for (srs::size_t i = 0; i < n; ++i) {
            b(i) = a(i);
        }
    }
}

// Copy n elements from a(j+i) to a(k+i).
template <class T>
void copy_n(srs::size_t n, srs::size_t j, srs::size_t k, Array<T, 1>& a)
{
    if (j > k) {
        for (srs::size_t i = 0; i < n; ++i) {
            a(k + i) = a(j + i);
        }
    }
    else if (j < k) {
        for (srs::size_t i = n; i > 0; --i) {
            a(k + i) = a(j + i);
        }
    }
}

template <class T>
void copy_n(srs::size_t n, srs::size_t j, srs::size_t k, Array_ref<T, 1>& a)
{
    if (j > k) {
        for (srs::size_t i = 0; i < n; ++i) {
            a(k + i) = a(j + i);
        }
    }
    else if (j < k) {
        for (srs::size_t i = n; i > 0; --i) {
            a(k + i) = a(j + i);
        }
    }
}

}  // namespace srs

#endif  // SRS_ARRAY_OPR_H
