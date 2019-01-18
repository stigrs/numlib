// Copyright (c) 2008-2010 Kent State University
// Copyright (c) 2011-2012 Texas A&M University
// Copyright (c) 2018 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#ifndef NUMLIB_DENSE_MATRIX_OPERATIONS_H
#define NUMLIB_DENSE_MATRIX_OPERATIONS_H

#ifdef USE_MKL
#include <mkl.h>
#else
#include <cblas.h>
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

#include <algorithm>
#include <random>
#include <iomanip>
#include <iostream>

namespace Numlib {

//------------------------------------------------------------------------------
//
// The following operations are defined for all Matrix types:

// Return matrix order.
template <typename M>
inline Enable_if<Matrix_type<M>(), std::size_t> rank(const M& m)
{
    return m.rank();
}

// Return matrix size.
template <typename M>
inline Enable_if<Matrix_type<M>(), Index> size(const M& m)
{
    return m.size();
}

// Return number of rows.
template <typename M>
inline Enable_if<Matrix_type<M>(), Index> rows(const M& m)
{
    static_assert(0 < M::order, "");
    return m.extent(0);
}

// Return number of columns.
template <typename M>
inline Enable_if<Matrix_type<M>(), Index> cols(const M& m)
{
    static_assert(1 < M::order, "");
    return m.extent(1);
}

// Return extent.
template <typename M>
inline Enable_if<Matrix_type<M>(), Index> extent(const M& m, Index dim)
{
    assert(dim < static_cast<Index>(m.order));
    return m.extent(dim);
}

// Create matrix of zeros.
template <typename M, typename... Args>
inline Enable_if<Matrix_type<M>(), M> zeros(Args... args)
{
    using value_type = typename M::value_type;

    assert(M::order == sizeof...(args));
    M res(args...);
    res = value_type{0};
    return res;
}

// Create matrix of ones.
template <typename M, typename... Args>
inline Enable_if<Matrix_type<M>(), M> ones(Args... args)
{
    using value_type = typename M::value_type;

    assert(M::order == sizeof...(args));
    M res(args...);
    res = value_type{1};
    return res;
}

// Create a random matrix from a normal distribution with zero mean and unit
// variance.
template <typename M, typename... Args>
inline Enable_if<Matrix_type<M>() && Real_type<typename M::value_type>(), M>
randn(Args... args)
{
    assert(M::order == sizeof...(args));
    M res(args...);

    std::random_device rd{};
    std::mt19937_64 gen{rd()};
    std::normal_distribution<> nd{};

    for (auto& x : res) {
        x = nd(gen);
    }
    return res;
}

// Create a random matrix from a uniform real distribution on the
// interval [0, 1).
template <typename M, typename... Args>
inline Enable_if<Matrix_type<M>() && Real_type<typename M::value_type>(), M>
randu(Args... args)
{
    assert(M::order == sizeof...(args));
    M res(args...);

    std::random_device rd{};
    std::mt19937_64 gen{rd()};

    std::uniform_real_distribution<> ur{};
    for (auto& x : res) {
        x = ur(gen);
    }
    return res;
}

// Create a random matrix from a uniform integer distribution on the
// interval [0, 1].
template <typename M, typename... Args>
inline Enable_if<Matrix_type<M>() && Integer_type<typename M::value_type>(), M>
randi(Args... args)
{
    assert(M::order == sizeof...(args));
    M res(args...);

    std::random_device rd{};
    std::mt19937_64 gen{rd()};

    std::uniform_int_distribution<> ui{};
    for (auto& x : res) {
        x = ui(gen);
    }
    return res;
}

//------------------------------------------------------------------------------
//
// Special methods for 2D matrices:

// Transpose.
template <typename T>
inline Matrix<T, 2> transpose(const Matrix<T, 2>& m)
{
    const Index n = m.rows();
    const Index p = m.cols();

    Matrix<T, 2> res(p, n);

#pragma omp parallel for
    for (Index i = 0; i < p; ++i) {
        for (Index j = 0; j < n; ++j) {
            res(i, j) = m.data()[i + j * p];
        }
    }
    return res;
}

//------------------------------------------------------------------------------
//
// Equality comparable:

// Two matrices compare equal when they have the same elements. Comparison
// of matrices decribed by different slices is undefined behavior.

template <typename M1, typename M2>
inline Enable_if<Matrix_type<M1>() && Matrix_type<M2>(), bool>
operator==(const M1& a, const M2& b)
{
    assert(same_extents(a, b));
    return std::equal(a.begin(), a.end(), b.begin());
}

template <typename M1, typename M2>
inline Enable_if<Matrix_type<M1>() && Matrix_type<M2>(), bool>
operator!=(const M1& a, const M2& b)
{
    return !(a == b);
}

template <typename M1, typename M2>
inline Enable_if<Matrix_type<M1>() && Matrix_type<M2>(), bool>
operator<(const M1& a, const M2& b)
{
    return std::lexicographical_compare(a.begin(), a.end(), b.begin(), b.end());
}

template <typename M1, typename M2>
inline Enable_if<Matrix_type<M1>() && Matrix_type<M2>(), bool>
operator>(const M1& a, const M2& b)
{
    return b < a;
}

template <typename M1, typename M2>
inline Enable_if<Matrix_type<M1>() && Matrix_type<M2>(), bool>
operator<=(const M1& a, const M2& b)
{
    return !(a > b);
}

template <typename M1, typename M2>
inline Enable_if<Matrix_type<M1>() && Matrix_type<M2>(), bool>
operator>=(const M1& a, const M2& b)
{
    return !(a < b);
}

//------------------------------------------------------------------------------
//
// Binary arithmetic operations:

// Matrix addition:

template <typename T, std::size_t N>
inline Matrix<T, N> operator+(const Matrix<T, N>& a, const Matrix<T, N>& b)
{
    assert(same_extents(a, b));

    Matrix<T, N> res = a;
    res += b;
    return res;
}

template <typename T, std::size_t N>
inline Matrix<T, N> operator+(const Matrix<const T, N>& a,
                              const Matrix<const T, N>& b)
{
    assert(same_extents(a, b));

    Matrix<T, N> res = a;
    res += b;
    return res;
}

template <typename T, std::size_t N>
inline Matrix<T, N> operator+(const Matrix_ref<T, N>& a,
                              const Matrix_ref<T, N>& b)
{
    assert(same_extents(a, b));

    Matrix<T, N> res = a;
    res += b;
    return res;
}

template <typename T, std::size_t N>
inline Matrix<T, N> operator+(const Matrix_ref<const T, N>& a,
                              const Matrix_ref<const T, N>& b)
{
    assert(same_extents(a, b));

    Matrix<T, N> res = a;
    res += b;
    return res;
}

template <typename T, std::size_t N>
inline Matrix<T, N> operator+(const Matrix<T, N>& a, const Matrix_ref<T, N>& b)
{
    assert(same_extents(a, b));

    Matrix<T, N> res = a;
    res += b;
    return res;
}

template <typename T, std::size_t N>
inline Matrix<T, N> operator+(const Matrix<const T, N>& a,
                              const Matrix_ref<const T, N>& b)
{
    assert(same_extents(a, b));

    Matrix<T, N> res = a;
    res += b;
    return res;
}

template <typename T, std::size_t N>
inline Matrix<T, N> operator+(const Matrix_ref<T, N>& a, const Matrix<T, N>& b)
{
    assert(same_extents(a, b));

    Matrix<T, N> res = a;
    res += b;
    return res;
}

template <typename T, std::size_t N>
inline Matrix<T, N> operator+(const Matrix_ref<const T, N>& a,
                              const Matrix<const T, N>& b)
{
    assert(same_extents(a, b));

    Matrix<T, N> res = a;
    res += b;
    return res;
}

template <typename T, std::size_t N>
inline Matrix<T, N> operator+(const Matrix_ref<T, N>& a,
                              const Matrix<const T, N>& b)
{
    assert(same_extents(a, b));

    Matrix<T, N> res = a;
    res += b;
    return res;
}

template <typename T, std::size_t N>
inline Matrix<T, N> operator+(const Matrix_ref<const T, N>& a,
                              const Matrix<T, N>& b)
{
    assert(same_extents(a, b));

    Matrix<T, N> res = a;
    res += b;
    return res;
}

template <typename T, std::size_t N>
inline Matrix<T, N> operator+(const Matrix<T, N>& a,
                              const Matrix_ref<const T, N>& b)
{
    assert(same_extents(a, b));

    Matrix<T, N> res = a;
    res += b;
    return res;
}

template <typename T, std::size_t N>
inline Matrix<T, N> operator+(const Matrix<const T, N>& a,
                              const Matrix_ref<T, N>& b)
{
    assert(same_extents(a, b));

    Matrix<T, N> res = a;
    res += b;
    return res;
}

// Matrix subtraction:

template <typename T, std::size_t N>
inline Matrix<T, N> operator-(const Matrix<T, N>& a, const Matrix<T, N>& b)
{
    assert(same_extents(a, b));

    Matrix<T, N> res = a;
    res -= b;
    return res;
}

template <typename T, std::size_t N>
inline Matrix<T, N> operator-(const Matrix<const T, N>& a,
                              const Matrix<const T, N>& b)
{
    assert(same_extents(a, b));

    Matrix<T, N> res = a;
    res -= b;
    return res;
}

template <typename T, std::size_t N>
inline Matrix<T, N> operator-(const Matrix_ref<T, N>& a,
                              const Matrix_ref<T, N>& b)
{
    assert(same_extents(a, b));

    Matrix<T, N> res = a;
    res -= b;
    return res;
}

template <typename T, std::size_t N>
inline Matrix<T, N> operator-(const Matrix_ref<const T, N>& a,
                              const Matrix_ref<const T, N>& b)
{
    assert(same_extents(a, b));

    Matrix<T, N> res = a;
    res -= b;
    return res;
}

template <typename T, std::size_t N>
inline Matrix<T, N> operator-(const Matrix<T, N>& a, const Matrix_ref<T, N>& b)
{
    assert(same_extents(a, b));

    Matrix<T, N> res = a;
    res -= b;
    return res;
}

template <typename T, std::size_t N>
inline Matrix<T, N> operator-(const Matrix<const T, N>& a,
                              const Matrix_ref<const T, N>& b)
{
    assert(same_extents(a, b));

    Matrix<T, N> res = a;
    res -= b;
    return res;
}

template <typename T, std::size_t N>
inline Matrix<T, N> operator-(const Matrix_ref<T, N>& a, const Matrix<T, N>& b)
{
    assert(same_extents(a, b));

    Matrix<T, N> res = a;
    res -= b;
    return res;
}

template <typename T, std::size_t N>
inline Matrix<T, N> operator-(const Matrix_ref<const T, N>& a,
                              const Matrix<const T, N>& b)
{
    assert(same_extents(a, b));

    Matrix<T, N> res = a;
    res -= b;
    return res;
}

template <typename T, std::size_t N>
inline Matrix<T, N> operator-(const Matrix_ref<T, N>& a,
                              const Matrix<const T, N>& b)
{
    assert(same_extents(a, b));

    Matrix<T, N> res = a;
    res -= b;
    return res;
}

template <typename T, std::size_t N>
inline Matrix<T, N> operator-(const Matrix_ref<const T, N>& a,
                              const Matrix<T, N>& b)
{
    assert(same_extents(a, b));

    Matrix<T, N> res = a;
    res -= b;
    return res;
}

template <typename T, std::size_t N>
inline Matrix<T, N> operator-(const Matrix<T, N>& a,
                              const Matrix_ref<const T, N>& b)
{
    assert(same_extents(a, b));

    Matrix<T, N> res = a;
    res -= b;
    return res;
}

template <typename T, std::size_t N>
inline Matrix<T, N> operator-(const Matrix<const T, N>& a,
                              const Matrix_ref<T, N>& b)
{
    assert(same_extents(a, b));

    Matrix<T, N> res = a;
    res -= b;
    return res;
}

// Scalar addition:

template <typename T, std::size_t N>
inline Matrix<T, N> operator+(const Matrix<T, N>& a, const T& scalar)
{
    Matrix<T, N> res = a;
    res += scalar;
    return res;
}

template <typename T, std::size_t N>
inline Matrix<T, N> operator+(const Matrix<const T, N>& a, const T& scalar)
{
    Matrix<T, N> res = a;
    res += scalar;
    return res;
}

template <typename T, std::size_t N>
inline Matrix<T, N> operator+(const Matrix_ref<T, N>& a, const T& scalar)
{
    Matrix<T, N> res = a;
    res += scalar;
    return res;
}

template <typename T, std::size_t N>
inline Matrix<T, N> operator+(const Matrix_ref<const T, N>& a, const T& scalar)
{
    Matrix<T, N> res = a;
    res += scalar;
    return res;
}

template <typename T, std::size_t N>
inline Matrix<T, N> operator+(const T& scalar, const Matrix<T, N>& a)
{
    Matrix<T, N> res = a;
    res += scalar;
    return res;
}

template <typename T, std::size_t N>
inline Matrix<T, N> operator+(const T& scalar, const Matrix<const T, N>& a)
{
    Matrix<T, N> res = a;
    res += scalar;
    return res;
}

template <typename T, std::size_t N>
inline Matrix<T, N> operator+(const T& scalar, const Matrix_ref<T, N>& a)
{
    Matrix<T, N> res = a;
    res += scalar;
    return res;
}

template <typename T, std::size_t N>
inline Matrix<T, N> operator+(const T& scalar, const Matrix_ref<const T, N>& a)
{
    Matrix<T, N> res = a;
    res += scalar;
    return res;
}

// Scalar subtraction:

template <typename T, std::size_t N>
inline Matrix<T, N> operator-(const Matrix<T, N>& a, const T& scalar)
{
    Matrix<T, N> res = a;
    res -= scalar;
    return res;
}

template <typename T, std::size_t N>
inline Matrix<T, N> operator-(const Matrix<const T, N>& a, const T& scalar)
{
    Matrix<T, N> res = a;
    res -= scalar;
    return res;
}

template <typename T, std::size_t N>
inline Matrix<T, N> operator-(const Matrix_ref<T, N>& a, const T& scalar)
{
    Matrix<T, N> res = a;
    res -= scalar;
    return res;
}

template <typename T, std::size_t N>
inline Matrix<T, N> operator-(const Matrix_ref<const T, N>& a, const T& scalar)
{
    Matrix<T, N> res = a;
    res -= scalar;
    return res;
}

// Scalar multiplication:

template <typename T, std::size_t N>
inline Matrix<T, N> operator*(const Matrix<T, N>& a, const T& scalar)
{
    Matrix<T, N> res = a;
    res *= scalar;
    return res;
}

template <typename T, std::size_t N>
inline Matrix<T, N> operator*(const Matrix<const T, N>& a, const T& scalar)
{
    Matrix<T, N> res = a;
    res *= scalar;
    return res;
}

template <typename T, std::size_t N>
inline Matrix<T, N> operator*(const Matrix_ref<T, N>& a, const T& scalar)
{
    Matrix<T, N> res = a;
    res *= scalar;
    return res;
}

template <typename T, std::size_t N>
inline Matrix<T, N> operator*(const Matrix_ref<const T, N>& a, const T& scalar)
{
    Matrix<T, N> res = a;
    res *= scalar;
    return res;
}

template <typename T, std::size_t N>
inline Matrix<T, N> operator*(const T& scalar, const Matrix<T, N>& a)
{
    Matrix<T, N> res = a;
    res *= scalar;
    return res;
}

template <typename T, std::size_t N>
inline Matrix<T, N> operator*(const T& scalar, const Matrix<const T, N>& a)
{
    Matrix<T, N> res = a;
    res *= scalar;
    return res;
}

template <typename T, std::size_t N>
inline Matrix<T, N> operator*(const T& scalar, const Matrix_ref<T, N>& a)
{
    Matrix<T, N> res = a;
    res *= scalar;
    return res;
}

template <typename T, std::size_t N>
inline Matrix<T, N> operator*(const T& scalar, const Matrix_ref<const T, N>& a)
{
    Matrix<T, N> res = a;
    res *= scalar;
    return res;
}

// Scalar division:

template <typename T, std::size_t N>
inline Matrix<T, N> operator/(const Matrix<T, N>& a, const T& scalar)
{
    Matrix<T, N> res = a;
    res /= scalar;
    return res;
}

template <typename T, std::size_t N>
inline Matrix<T, N> operator/(const Matrix<const T, N>& a, const T& scalar)
{
    Matrix<T, N> res = a;
    res /= scalar;
    return res;
}

template <typename T, std::size_t N>
inline Matrix<T, N> operator/(const Matrix_ref<T, N>& a, const T& scalar)
{
    Matrix<T, N> res = a;
    res /= scalar;
    return res;
}

template <typename T, std::size_t N>
inline Matrix<T, N> operator/(const Matrix_ref<const T, N>& a, const T& scalar)
{
    Matrix<T, N> res = a;
    res /= scalar;
    return res;
}

// Scalar modulus:

template <typename T, std::size_t N>
inline Enable_if<Integer_type<T>(), Matrix<T, N>>
operator%(const Matrix<T, N>& a, const T& scalar)
{
    Matrix<T, N> res = a;
    res %= scalar;
    return res;
}

template <typename T, std::size_t N>
inline Enable_if<Integer_type<T>(), Matrix<T, N>>
operator%(const Matrix<const T, N>& a, const T& scalar)
{
    Matrix<T, N> res = a;
    res %= scalar;
    return res;
}

template <typename T, std::size_t N>
inline Enable_if<Integer_type<T>(), Matrix<T, N>>
operator%(const Matrix_ref<T, N>& a, const T& scalar)
{
    Matrix<T, N> res = a;
    res %= scalar;
    return res;
}

template <typename T, std::size_t N>
inline Enable_if<Integer_type<T>(), Matrix<T, N>>
operator%(const Matrix_ref<const T, N>& a, const T& scalar)
{
    Matrix<T, N> res = a;
    res %= scalar;
    return res;
}

//------------------------------------------------------------------------------
//
// Matrix-matrix multiplication:

// Multiplication of N x M by M x P matrix.
template <typename M1, typename M2, typename M3>
Enable_if<Matrix_type<M1>() && Matrix_type<M2>() && Matrix_type<M3>(), void>
mm_mul(const M1& a, const M2& b, M3& res)
{
    static_assert(M1::order == 2, "bad order for matrix-matrix multiplication");
    static_assert(M2::order == 2, "bad order for matrix-matrix multiplication");
    static_assert(M3::order == 2, "bad order for matrix-matrix multiplication");

    using value_type = typename M1::value_type;

    const Index n = a.extent(0);
    const Index m = a.extent(1);
    const Index p = b.extent(1);
    assert(m == b.extent(0));

    res.resize(n, p);

#pragma omp parallel for shared(res, a, b)
    for (Index i = 0; i < n; ++i) {
        for (Index j = 0; j < p; ++j) {
            res(i, j) = value_type{0};
            for (Index k = 0; k < m; ++k) {
                res(i, j) += a(i, k) * b(k, j);
            }
        }
    }
}

// Use BLAS for double matrices.
inline void mm_mul(const Matrix<double, 2>& a,
                   const Matrix<double, 2>& b,
                   Matrix<double, 2>& res)
{
    constexpr double alpha = 1.0;
    constexpr double beta = 0.0;

    const BLAS_INT m = narrow_cast<BLAS_INT>(a.rows());
    const BLAS_INT n = narrow_cast<BLAS_INT>(b.cols());
    const BLAS_INT k = narrow_cast<BLAS_INT>(a.cols());

    const BLAS_INT lda = k;
    const BLAS_INT ldb = n;
    const BLAS_INT ldc = n;

    res.resize(m, n);

    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k, alpha,
                a.data(), lda, b.data(), ldb, beta, res.data(), ldc);
}

template <typename T>
inline Matrix<T, 2> operator*(const Matrix<T, 2>& a, const Matrix<T, 2>& b)
{
    Matrix<T, 2> res;
    mm_mul(a, b, res);
    return res;
}

template <typename T>
inline Matrix<T, 2> operator*(const Matrix<const T, 2>& a,
                              const Matrix<const T, 2>& b)
{
    Matrix<T, 2> res;
    mm_mul(a, b, res);
    return res;
}

template <typename T>
inline Matrix<T, 2> operator*(const Matrix_ref<T, 2>& a,
                              const Matrix_ref<T, 2>& b)
{
    Matrix<T, 2> res;
    mm_mul(a, b, res);
    return res;
}

template <typename T>
inline Matrix<T, 2> operator*(const Matrix_ref<const T, 2>& a,
                              const Matrix_ref<const T, 2>& b)
{
    Matrix<T, 2> res;
    mm_mul(a, b, res);
    return res;
}

template <typename T>
inline Matrix<T, 2> operator*(const Matrix<T, 2>& a, const Matrix_ref<T, 2>& b)
{
    Matrix<T, 2> res;
    mm_mul(a, b, res);
    return res;
}

template <typename T>
inline Matrix<T, 2> operator*(const Matrix<const T, 2>& a,
                              const Matrix_ref<const T, 2>& b)
{
    Matrix<T, 2> res;
    mm_mul(a, b, res);
    return res;
}

template <typename T>
inline Matrix<T, 2> operator*(const Matrix_ref<T, 2>& a, const Matrix<T, 2>& b)
{
    Matrix<T, 2> res;
    mm_mul(a, b, res);
    return res;
}

template <typename T>
inline Matrix<T, 2> operator*(const Matrix_ref<const T, 2>& a,
                              const Matrix<const T, 2>& b)
{
    Matrix<T, 2> res;
    mm_mul(a, b, res);
    return res;
}

//------------------------------------------------------------------------------
//
// Matrix-vector multiplication:

template <typename M1, typename M2, typename M3>
Enable_if<Matrix_type<M1>() && Matrix_type<M2>() && Matrix_type<M3>(), void>
mv_mul(const M1& a, const M2& x, M3& y)
{
    static_assert(M1::order == 2, "bad order for matrix-vector multiplication");
    static_assert(M2::order == 1, "bad order for matrix-vector multiplication");
    static_assert(M3::order == 1, "bad order for matrix-vector multiplication");

    assert(x.size() == a.cols());

    using value_type = typename M1::value_type;

    y.resize(a.rows());

#pragma omp parallel for shared(y, a, x)
    for (Index i = 0; i < a.rows(); ++i) {
        y(i) = value_type{0};
        for (Index j = 0; j < a.cols(); ++j) {
            y(i) += a(i, j) * x(j);
        }
    }
}

// Use BLAS for double matrices and vectors.
inline void mv_mul(const Matrix<double, 2>& a,
                   const Matrix<double, 1>& x,
                   Matrix<double, 1>& y)
{
    constexpr double alpha = 1.0;
    constexpr double beta = 0.0;

    assert(x.size() == a.cols());

    const BLAS_INT m = narrow_cast<BLAS_INT>(a.rows());
    const BLAS_INT n = narrow_cast<BLAS_INT>(a.cols());

    y.resize(m);

    const BLAS_INT lda = n;
    const BLAS_INT incx = narrow_cast<BLAS_INT>(x.descriptor().strides[0]);
    const BLAS_INT incy = narrow_cast<BLAS_INT>(y.descriptor().strides[0]);

    cblas_dgemv(CblasRowMajor, CblasNoTrans, m, n, alpha, a.data(), lda,
                x.data(), incx, beta, y.data(), incy);
}

template <typename T>
inline Matrix<T, 1> operator*(const Matrix<T, 2>& a, const Matrix<T, 1>& x)
{
    Matrix<T, 1> res;
    mv_mul(a, x, res);
    return res;
}

template <typename T>
inline Matrix<T, 1> operator*(const Matrix<const T, 2>& a,
                              const Matrix<const T, 1>& x)
{
    Matrix<T, 1> res;
    mv_mul(a, x, res);
    return res;
}

template <typename T>
inline Matrix<T, 1> operator*(const Matrix_ref<T, 2>& a,
                              const Matrix_ref<T, 1>& x)
{
    Matrix<T, 1> res;
    mv_mul(a, x, res);
    return res;
}

template <typename T>
inline Matrix<T, 1> operator*(const Matrix_ref<const T, 2>& a,
                              const Matrix_ref<const T, 1>& x)
{
    Matrix<T, 1> res;
    mv_mul(a, x, res);
    return res;
}

template <typename T>
inline Matrix<T, 1> operator*(const Matrix<T, 2>& a, const Matrix_ref<T, 1>& x)
{
    Matrix<T, 1> res;
    mv_mul(a, x, res);
    return res;
}

template <typename T>
inline Matrix<T, 1> operator*(const Matrix<const T, 2>& a,
                              const Matrix_ref<const T, 1>& x)
{
    Matrix<T, 1> res;
    mv_mul(a, x, res);
    return res;
}

template <typename T>
inline Matrix<T, 1> operator*(const Matrix_ref<T, 2>& a, const Matrix<T, 1>& x)
{
    Matrix<T, 1> res;
    mv_mul(a, x, res);
    return res;
}

template <typename T>
inline Matrix<T, 1> operator*(const Matrix_ref<const T, 2>& a,
                              const Matrix<const T, 1>& x)
{
    Matrix<T, 1> res;
    mv_mul(a, x, res);
    return res;
}

//------------------------------------------------------------------------------
//
// Hadamard product:
//
// The hadamard product can easily be generalized to N-dimensional matrices
// since the operation is performed element-wise. The operands only need
// to be of the same shape.

template <typename M1, typename M2, typename M3>
Enable_if<Matrix_type<M1>() && Matrix_type<M2>() && Matrix_type<M3>(), void>
hadamard_product(const M1& a, const M2& b, M3& res)
{
    static_assert((M1::order == M2::order) && (M2::order == M3::order),
                  "bad matrix order for hadamard product");
    assert(a.shape() == b.shape());
    res.resize(a.shape());

    using Mul = std::multiplies<Value_type<M1>>;
    std::transform(a.begin(), a.end(), b.begin(), res.begin(), Mul{});
}

//------------------------------------------------------------------------------
//
// I/O operators for 1D and 2D matrices:

template <typename T>
std::ostream& operator<<(std::ostream& to, const Matrix<T, 1>& a)
{
    to << a.size() << '\n' << "[ ";
    for (Index i = 0; i < a.size(); ++i) {
        to << std::setw(9) << a(i) << " ";
        if (!((i + 1) % 7) && (i != (a.size() - 1))) {
            to << "\n  ";
        }
    }
    to << ']';
    return to;
}

template <typename T>
std::ostream& operator<<(std::ostream& to, const Matrix_ref<T, 1>& a)
{
    to << a.size() << '\n' << "[ ";
    for (Index i = 0; i < a.size(); ++i) {
        to << std::setw(9) << a(i) << " ";
        if (!((i + 1) % 7) && (i != (a.size() - 1))) {
            to << "\n  ";
        }
    }
    to << ']';
    return to;
}

// The vector must be entered as: n [...].
template <typename T>
std::istream& operator>>(std::istream& from, Matrix<T, 1>& a)
{
    Index n;
    from >> n;
    a.resize(n);

    char ch;
    from >> ch; // [
    for (Index i = 0; i < n; ++i) {
        from >> a(i);
    }
    from >> ch; // ]
    return from;
}

template <typename T>
std::ostream& operator<<(std::ostream& to, const Matrix<T, 2>& a)
{
    to << a.rows() << " x " << a.cols() << "\n[";
    for (Index i = 0; i < a.rows(); ++i) {
        for (Index j = 0; j < a.cols(); ++j) {
            to << std::setw(9) << a(i, j) << " ";
        }
        if (i != (a.rows() - 1)) {
            to << "\n ";
        }
    }
    to << "]\n";
    return to;
}

template <typename T>
std::ostream& operator<<(std::ostream& to, const Matrix_ref<T, 2>& a)
{
    to << a.rows() << " x " << a.cols() << "\n[";
    for (Index i = 0; i < a.rows(); ++i) {
        for (Index j = 0; j < a.cols(); ++j) {
            to << std::setw(9) << a(i, j) << " ";
        }
        if (i != (a.rows() - 1)) {
            to << "\n ";
        }
    }
    to << "]\n";
    return to;
}

// The matrix must be entered as: m x n [...].
template <typename T>
std::istream& operator>>(std::istream& from, Matrix<T, 2>& a)
{
    Index m;
    Index n;
    char ch;

    from >> m >> ch >> n; // m x n
    a.resize(m, n);

    from >> ch; // [
    for (Index i = 0; i < m; ++i) {
        for (Index j = 0; j < n; ++j) {
            from >> a(i, j);
        }
    }
    from >> ch; // ]
    return from;
}

// Output stream operator for 3D matrices:

template <typename T>
std::ostream& operator<<(std::ostream& to, const Matrix<T, 3>& a)
{
    to << a.extent(0) << " x " << a.extent(1) << " x " << a.extent(2) << "\n[";
    for (Index k = 0; k < a.extent(2); ++k) {
        for (Index i = 0; i < a.extent(0); ++i) {
            for (Index j = 0; j < a.extent(1); ++j) {
                to << std::setw(9) << a(i, j, k) << " ";
            }
            if (i != (a.extent(0) - 1)) {
                to << "\n ";
            }
        }
        if (k != (a.extent(2) - 1)) {
            to << "\n\n ";
        }
    }
    to << "]\n";
    return to;
}

template <typename T>
std::ostream& operator<<(std::ostream& to, const Matrix_ref<T, 3>& a)
{
    to << a.extent(0) << " x " << a.extent(1) << " x " << a.extent(2) << "\n[";
    for (Index k = 0; k < a.extent(2); ++k) {
        for (Index i = 0; i < a.extent(0); ++i) {
            for (Index j = 0; j < a.extent(1); ++j) {
                to << std::setw(9) << a(i, j, k) << " ";
            }
            if (i != (a.extent(0) - 1)) {
                to << "\n ";
            }
        }
        if (k != (a.extent(2) - 1)) {
            to << "\n\n ";
        }
    }
    to << "]\n";
    return to;
}

} // namespace Numlib

#endif // NUMLIB_DENSE_MATRIX_OPERATIONS_H
