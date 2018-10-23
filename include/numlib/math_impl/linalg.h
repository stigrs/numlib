// Copyright (c) 2018 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#ifndef NUMLIB_MATH_LINALG_H
#define NUMLIB_MATH_LINALG_H

#ifdef USE_MKL
#include <mkl.h>
#else
#include <cblas.h>
#include <lapacke.h>
#endif

#include <numlib/traits.h>
#include <numlib/matrix.h>
#include <algorithm>
#include <numeric>
#include <functional>
#include <cmath>
#include <complex>

namespace Numlib {

//------------------------------------------------------------------------------
//
// Create special vectors and matrices:

// Create linearly spaced vector.
Vec<double> linspace(double x1, double x2, Index n);

// Identity matrix.
template <typename T = double>
inline Mat<T> identity(Index n)
{
    Mat<T> res = zeros<Mat<T>>(n, n);
    res.diag() = T{1};
    return res;
}

// Hilbert matrix.
template <typename T = double>
inline Enable_if<Real_type<T>(), Mat<T>> hilbert(Index n)
{
    Mat<T> res(n, n);
    for (Index i = 0; i < n; ++i) {
        for (Index j = 0; j < n; ++j) {
            res(i, j) = T{1} / static_cast<T>(i + j + 1);
        }
    }
    return res;
}

//------------------------------------------------------------------------------
//
// Find max, min, sum, and product of elements:

template <typename M>
inline Enable_if<Matrix_type<M>(), typename M::value_type> max(const M& vec)
{
    static_assert(M::order == 1, "bad rank for max(vec)");
    return *std::max_element(vec.begin(), vec.end());
}

template <typename M>
inline Enable_if<Matrix_type<M>(), Vec<typename M::value_type>>
max(const M& mat, Index dim)
{
    static_assert(M::order == 2, "bad rank for max(mat)");
    assert(dim >= 0 && dim < static_cast<Index>(M::order));

    Vec<Value_type<M>> result(mat.extent(dim));
    if (dim == 0) { // row
        for (Index i = 0; i < mat.rows(); ++i) {
            result(i) = max(mat.row(i));
        }
    }
    else { // column
        for (Index i = 0; i < mat.cols(); ++i) {
            result(i) = max(mat.column(i));
        }
    }
    return result;
}

template <typename M>
inline Enable_if<Matrix_type<M>(), typename M::value_type> min(const M& vec)
{
    static_assert(M::order == 1, "bad rank for min(vec)");
    return *std::min_element(vec.begin(), vec.end());
}

template <typename M>
inline Enable_if<Matrix_type<M>(), Vec<typename M::value_type>>
min(const M& mat, Index dim)
{
    static_assert(M::order == 2, "bad rank for min(mat)");
    assert(dim >= 0 && dim < static_cast<Index>(M::order));

    Vec<Value_type<M>> result(mat.extent(dim));
    if (dim == 0) { // row
        for (Index i = 0; i < mat.rows(); ++i) {
            result(i) = min(mat.row(i));
        }
    }
    else { // column
        for (Index i = 0; i < mat.cols(); ++i) {
            result(i) = min(mat.column(i));
        }
    }
    return result;
}

template <typename M>
inline Enable_if<Matrix_type<M>(), typename M::value_type> sum(const M& vec)
{
    static_assert(M::order == 1, "bad rank for sum(vec)");
    constexpr auto zero = Value_type<M>{0};
    return std::accumulate(vec.begin(), vec.end(), zero);
}

template <typename M>
inline Enable_if<Matrix_type<M>(), Vec<typename M::value_type>>
sum(const M& mat, Index dim)
{
    static_assert(M::order == 2, "bad rank for sum(mat)");
    assert(dim >= 0 && dim < static_cast<Index>(M::order));

    Vec<Value_type<M>> result(mat.extent(dim));
    if (dim == 0) { // row
        for (Index i = 0; i < mat.rows(); ++i) {
            result(i) = sum(mat.row(i));
        }
    }
    else { // column
        for (Index i = 0; i < mat.cols(); ++i) {
            result(i) = sum(mat.column(i));
        }
    }
    return result;
}

template <typename M>
inline Enable_if<Matrix_type<M>(), typename M::value_type> prod(const M& vec)
{
    static_assert(M::order == 1, "bad rank for prod(vec)");

    using T = typename M::value_type;
    constexpr auto one = T{1};
    return std::accumulate(vec.begin(), vec.end(), one, std::multiplies<T>());
}

template <typename M>
inline Enable_if<Matrix_type<M>(), Vec<typename M::value_type>>
prod(const M& mat, Index dim)
{
    static_assert(M::order == 2, "bad rank for prod(mat)");
    assert(dim >= 0 && dim < static_cast<Index>(M::order));

    Vec<Value_type<M>> result(mat.extent(dim));
    if (dim == 0) { // row
        for (Index i = 0; i < mat.rows(); ++i) {
            result(i) = prod(mat.row(i));
        }
    }
    else { // column
        for (Index i = 0; i < mat.cols(); ++i) {
            result(i) = prod(mat.column(i));
        }
    }
    return result;
}

//------------------------------------------------------------------------------
//
// Matrix and vector products:

// Dot product of dense vectors.
template <typename M1, typename M2>
inline Enable_if<Matrix_type<M1>() && Matrix_type<M2>(),
                 typename M1::value_type>
dot(const M1& x, const M2& y)
{
    static_assert(M1::order == 1, "bad rank for dot product");
    static_assert(M2::order == 1, "bad rank for dot product");
    assert(same_extents(x, y));

    constexpr auto zero = Value_type<M1>{0};
    return std::inner_product(x.begin(), x.end(), y.begin(), zero);
}

// Dot product of sparse vectors.
template <typename T>
inline T dot(const Sparse_vector<T>& x, const Sparse_vector<T>& y)
{
    assert(x.size() == y.size());

    T result{0};
    for (Index i = 0; i < x.size(); ++i) { // inefficient for large vectors
        result += x(i) * y(i);
    }
    return result;
}

// Dot product of a sparse and a dense vector.
template <typename T>
inline T dot(const Sparse_vector<T>& x, const Vec<T>& y)
{
    T result{0};

    Index i = 0;
    for (const auto& v : x) {
        result += v * y(x.loc(i));
        ++i;
    }
    return result;
}

// Dot product of a dense and a sparse vector.
template <typename T>
inline T dot(const Vec<T>& y, const Sparse_vector<T>& x)
{
    return dot(x, y);
}

// Cross product.
template <typename T>
inline void cross(const Vec<T>& x, const Vec<T>& y, Vec<T>& res)
{
    assert(x.size() == 3 && x.size() == y.size());
    res.resize(3);
    res(0) = x(1) * y(2) - x(2) * y(1);
    res(1) = x(2) * y(0) - x(0) * y(2);
    res(2) = x(0) * y(1) - x(1) * y(0);
}

// Cross product.
template <typename M1, typename M2>
inline Enable_if<Matrix_type<M1>() && Matrix_type<M2>(),
                 Vec<typename M1::value_type>>
cross(const M1& x, const M2& y)
{
    static_assert(Same<Value_type<M1>(), Value_type<M2>()>(),
                  "cross: different value types");

    Vec<typename M1::value_type> res;
    cross(x, y, res);
    return res;
}

// Compute vector-scalar product and add the result to a vector.
template <typename T>
inline void axpy(const T& a, const Vec<T>& x, Vec<T>& y)
{
    assert(same_extents(x, y));
    for (Index i = 0; i < x.size(); ++i) {
        y(i) = a * x(i) + y(i);
    }
}

// Matrix-matrix multiplication.
inline void matmul(const Mat<double>& a, const Mat<double>& b, Mat<double>& res)
{
    mm_mul(a, b, res);
}

// Matrix-vector multiplication.
inline void matmul(const Mat<double>& a, const Vec<double>& x, Vec<double>& y)
{
    mv_mul(a, x, y);
}

// Kronecker product.
template <typename T>
void kron(const Mat<T>& a, const Mat<T>& b, Mat<T>& res)
{
    const Index m = a.rows();
    const Index n = a.cols();
    const Index p = b.rows();
    const Index q = b.cols();

    res.resize(m * p, n * q);

    for (Index i = 0; i < m; ++i) {
        for (Index j = 0; j < n; ++j) {
            Index i0 = i * p;
            Index j0 = j * q;
            for (Index k = 0; k < p; ++k) {
                for (Index l = 0; l < q; ++l) {
                    res(i0 + k, j0 + l) = a(i, j) * b(k, l);
                }
            }
        }
    }
}

//------------------------------------------------------------------------------
//
// Matrix decompositions:

// LU factorization.
inline void lu(Mat<double>& a, Vec<BLAS_INT>& ipiv)
{
    const BLAS_INT m = narrow_cast<BLAS_INT>(a.rows());
    const BLAS_INT n = narrow_cast<BLAS_INT>(a.cols());
    const BLAS_INT lda = n;

    ipiv.resize(std::min(m, n));

    BLAS_INT info =
        LAPACKE_dgetrf(LAPACK_ROW_MAJOR, m, n, a.data(), lda, ipiv.data());
    if (info < 0) {
        throw Math_error("dgetrf: illegal input parameter");
    }
    if (info > 0) {
        throw Math_error("dgetrf: U matrix is singular");
    }
}

// QR factorization.
inline void qr(const Mat<double>& a, Mat<double>& q, Mat<double>& r)
{
    // Compute QR factorization:

    const BLAS_INT m = narrow_cast<BLAS_INT>(a.rows());
    const BLAS_INT n = narrow_cast<BLAS_INT>(a.cols());
    const BLAS_INT lda = n;

    Vec<double> tau(std::min(m, n));

    q.resize(m, n);
    q = a;

    BLAS_INT info =
        LAPACKE_dgeqrf(LAPACK_ROW_MAJOR, m, n, q.data(), lda, tau.data());
    if (info != 0) {
        throw Math_error("dgeqrf failed");
    }

    // Compute Q:

    info = LAPACKE_dorgqr(LAPACK_ROW_MAJOR, m, n, n, q.data(), lda, tau.data());
    if (info != 0) {
        throw Math_error("dorgqr failed");
    }

    // Compute R:

    r.resize(m, n);
    mm_mul(transpose(q), a, r);
}

// Singular value decomposition.
inline void svd(Mat<double>& a, Vec<double>& s, Mat<double>& u, Mat<double>& vt)
{
    BLAS_INT m = narrow_cast<BLAS_INT>(a.rows());
    BLAS_INT n = narrow_cast<BLAS_INT>(a.cols());
    BLAS_INT lda = n;
    BLAS_INT ldu = m;
    BLAS_INT ldvt = n;

    s.resize(std::min(m, n));
    u.resize(m, ldu);
    vt.resize(n, ldvt);

    Vec<double> superb(std::min(m, n) - 1);

    BLAS_INT info =
        LAPACKE_dgesvd(LAPACK_ROW_MAJOR, 'A', 'A', m, n, a.data(), lda,
                       s.data(), u.data(), ldu, vt.data(), ldvt, superb.data());
    if (info != 0) {
        throw Math_error("dgesvd failed");
    }
}

// Schmidt orthogonalization of n orbitals in a.
void schmidt(Mat<double>& a, Index n);

//------------------------------------------------------------------------------
//
// Matrix eigenvalues:

// Compute eigenvalues and eigenvectors of a real symmetric matrix.
inline void eigs(Mat<double>& a, Vec<double>& w)
{
    assert(a.rows() == a.cols());

    const BLAS_INT n = narrow_cast<BLAS_INT>(a.rows());
    w.resize(n);

    BLAS_INT info =
        LAPACKE_dsyevd(LAPACK_ROW_MAJOR, 'V', 'U', n, a.data(), n, w.data());
    if (info != 0) {
        throw Math_error("dsyevd failed");
    }
}

// Compute eigenvalues and eigenvectors of a real non-symmetric matrix.
void eig(Mat<double>& a,
         Mat<std::complex<double>>& evec,
         Vec<std::complex<double>>& eval);

// Compute eigenvalues and eigenvectors of a real symmetric band matrix.
inline void eigs(Band_mat<double>& ab, Mat<double>& evec, Vec<double>& eval)
{
    assert(ab.rows() == ab.cols());
    assert(ab.upper() == ab.lower());

    evec.resize(ab.rows(), ab.cols());
    eval.resize(ab.cols());

    const BLAS_INT n = narrow_cast<BLAS_INT>(ab.cols());
    const BLAS_INT kd = narrow_cast<BLAS_INT>(ab.upper());
    const BLAS_INT ldab = narrow_cast<BLAS_INT>(ab.leading_dim());
    const BLAS_INT ldz = narrow_cast<BLAS_INT>(ab.cols());

    BLAS_INT info = LAPACKE_dsbev(LAPACK_COL_MAJOR, 'V', 'U', n, kd, ab.data(),
                                  ldab, eval.data(), evec.data(), ldz);
    if (info != 0) {
        throw Math_error("dsbev failed");
    }
}

// Compute eigenvalues and eigenvectors in the interval [emin, emax] for
// a real band matrix.
//
// Note:
// - Only Intel MKL is supported.
//
#ifdef USE_MKL
void eig(double emin,
         double emax,
         const Band_mat<double>& ab,
         Mat<double>& evec,
         Vec<double>& eval);
#endif

// Compute eigenvalues and eigenvectors of a real symmetric matrix held in
// packed storage.
//
// Note:
// - Only Intel MKL is supported on Windows since OpenBLAS v0.2.14.1 gives 
//   wrong results.
//
#if defined(__unix__) || defined(__APPLE__) || defined(USE_MKL)
template <Uplo_scheme Uplo>
void eigs(Symm_mat<double, Uplo>& ap, Mat<double>& evec, Vec<double>& eval)
{
    assert(ap.size() >= eval.size() * (eval.size() + 1) / 2);

    const BLAS_INT n = narrow_cast<BLAS_INT>(eval.size());
    const BLAS_INT ldz = n;
    char uplo = ap.uplo_scheme();

    evec.resize(n, n);

    BLAS_INT info = LAPACKE_dspevd(LAPACK_ROW_MAJOR, 'V', uplo, n, ap.data(),
                                   eval.data(), evec.data(), ldz);
    if (info != 0) {
        throw Math_error("dspevd failed");
    }
}
#endif

// Compute eigenvalues and eigenvectors in the interval [emin, emax] for
// a real sparse matrix.
#ifdef USE_MKL
void eig(double emin,
         double emax,
         const Sp_mat<double>& a,
         Mat<double>& evec,
         Vec<double>& eval);
#endif

//------------------------------------------------------------------------------
//
// Norms and other numbers:

// Euclidean norm of dense vector.
template <typename M>
inline Enable_if<Matrix_type<M>() && Real_type<Value_type<M>>(),
                 typename M::value_type>
norm(const M& vec)
{
    static_assert(M::order == 1, "norm: bad matrix rank");

    using T = typename M::value_type;

    T result = T{0};
    if (!vec.empty()) {
        for (const auto& x : vec) {
            result += x * x;
        }
        result = std::sqrt(result);
    }
    return result;
}

// Euclidean norm of sparse vector.
template <typename T>
inline T norm(const Sparse_vector<T>& vec)
{
    T result{0};
    if (!vec.empty()) {
        for (const auto& x : vec) {
            result += x * x;
        }
        result = std::sqrt(result);
    }
    return result;
}

// Normalize dense vector.
template <typename M>
inline Enable_if<Matrix_type<M>() && Real_type<Value_type<M>>(), M>
normalize(const M& vec)
{
    static_assert(M::order == 1, "normalize: bad matrix rank");
    constexpr auto zero = Value_type<M>{0};

    M result(vec);
    auto n = norm(vec);
    if (n > zero) {
        result /= n;
    }
    return result;
}

// Normalize sparse vector.
template <typename T>
inline Sparse_vector<T> normalize(const Sparse_vector<T>& vec)
{
    constexpr auto zero = T{0};

    Sparse_vector<T> result(vec);
    auto n = norm(vec);
    if (n > zero) {
        result /= n;
    }
    return result;
}

// Matrix norm of a general rectangular matrix:
//
// Types of matrix norms:
// - M, m:       largest absolute value of the matrix
// - 1, O, o:    1-norm of the matrix (maximum column sum)
// - I, i:       infinity norm of the matrix (maximum row sum)
// - F, f, E, e: Frobenius norm of the matrix (square root of sum of squares)
//
inline double norm(const Mat<double>& a, char norm)
{
    assert(norm == 'M' || norm == 'm' || norm == '1' || norm == 'O' ||
           norm == 'o' || norm == 'I' || norm == 'i' || norm == 'F' ||
           norm == 'f' || norm == 'E' || norm == 'e');

    BLAS_INT m = narrow_cast<BLAS_INT>(a.rows());
    BLAS_INT n = narrow_cast<BLAS_INT>(a.cols());
    BLAS_INT lda = n;

    return LAPACKE_dlange(LAPACK_ROW_MAJOR, norm, m, n, a.data(), lda);
}

// Trace of a square matrix.
template <typename M>
inline Enable_if<Matrix_type<M>(), typename M::value_type> trace(const M& mat)
{
    static_assert(M::order == 2, "trace: bad matrix rank");
    assert(mat.rows() == mat.cols());

    constexpr auto zero = Value_type<M>{0};

    const auto d = mat.diag();
    return std::accumulate(d.begin(), d.end(), zero);
}

// Determinant of square matrix.
double det(const Mat<double>& a);

// Reciprocal condition number of a matrix.
inline double rcond(const Mat<double>& a)
{
    char nrm = '1';
    double anorm = norm(a, nrm);

    Mat<double> tmp(a);
    Vec<BLAS_INT> ipiv;
    lu(tmp, ipiv);

    BLAS_INT m = narrow_cast<BLAS_INT>(a.rows());
    BLAS_INT lda = narrow_cast<BLAS_INT>(a.cols());

    double res;
    BLAS_INT info =
        LAPACKE_dgecon(LAPACK_ROW_MAJOR, nrm, m, tmp.data(), lda, anorm, &res);
    if (info != 0) {
        throw Math_error("dgecon failed");
    }
    return res;
}

// Condition number of a matrix.
inline double cond(const Mat<double>& a) { return 1.0 / rcond(a); }

//------------------------------------------------------------------------------
//
// Solving equations and inverting matrices:

// Solve linear system of equations.
inline void linsolve(Mat<double>& a, Mat<double>& b)
{
    assert(a.rows() == a.cols());
    assert(b.rows() == a.cols());

    const BLAS_INT n = narrow_cast<BLAS_INT>(a.cols());
    const BLAS_INT nrhs = narrow_cast<BLAS_INT>(b.cols());
    const BLAS_INT lda = n;
    const BLAS_INT ldb = nrhs;

    Vec<BLAS_INT> ipiv(n);

    BLAS_INT info = LAPACKE_dgesv(LAPACK_ROW_MAJOR, n, nrhs, a.data(), lda,
                                  ipiv.data(), b.data(), ldb);
    if (info != 0) {
        throw Math_error("dgesv: factor U is singular");
    }
}

// Solve linear system of equations for a real, nonsymmetric sparse matrix.
#ifdef USE_MKL
void linsolve(const Sp_mat<double>& a, Mat<double>& b, Mat<double>& x);
#endif

// Matrix inversion.
inline void inv(Mat<double>& a)
{
    assert(a.rows() == a.cols());

    if (det(a) == 0.0) {
        throw Math_error("inv: matrix not invertible");
    }
    const BLAS_INT n = narrow_cast<BLAS_INT>(a.rows());
    const BLAS_INT lda = n;

    Vec<BLAS_INT> ipiv(n);
    lu(a, ipiv); // perform LU factorization

    BLAS_INT info =
        LAPACKE_dgetri(LAPACK_ROW_MAJOR, n, a.data(), lda, ipiv.data());
    if (info != 0) {
        throw Math_error("dgetri: matrix inversion failed");
    }
}

// Compute the minimum norm-solution to a real linear least squares problem.
inline void lstsq(Mat<double>& a, Mat<double>& b)
{
    BLAS_INT m = narrow_cast<BLAS_INT>(a.rows());
    BLAS_INT n = narrow_cast<BLAS_INT>(a.cols());
    BLAS_INT nrhs = narrow_cast<BLAS_INT>(b.cols());
    BLAS_INT lda = n;
    BLAS_INT ldb = nrhs;
    BLAS_INT rank;

    double rcond = -1.0;           // use machine epsilon
    Vec<double> s(std::min(m, n)); // singular values of a

    BLAS_INT info = LAPACKE_dgelsd(LAPACK_ROW_MAJOR, m, n, nrhs, a.data(), lda,
                                   b.data(), ldb, s.data(), rcond, &rank);
    if (info != 0) {
        throw Math_error("dgelsd failed");
    }
}

} // namespace Numlib

#endif // NUMLIB_MATH_LINALG_H
