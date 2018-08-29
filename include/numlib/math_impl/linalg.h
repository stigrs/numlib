// Copyright (c) 2018 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#ifndef NUMLIB_MATH_LINALG_H
#define NUMLIB_MATH_LINALG_H

#include <numlib/matrix.h>
#include <numlib/traits.h>
#include <algorithm>
#include <numeric>
#include <functional>
#include <cmath>

#ifdef __APPLE__
#include <Accelerate/Accelerate.h>
#else
#include <cblas.h>
#endif

namespace Numlib {

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
max(const M& mat, std::size_t dim)
{
    static_assert(M::order == 2, "bad rank for max(mat)");
    assert(dim >= 0 && dim < M::order);

    Vec<Value_type<M>> result(mat.extent(dim));
    if (dim == 0) { // row
        for (std::size_t i = 0; i < mat.rows(); ++i) {
            result(i) = max(mat.row(i));
        }
    }
    else { // column
        for (std::size_t i = 0; i < mat.rows(); ++i) {
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
min(const M& mat, std::size_t dim)
{
    static_assert(M::order == 2, "bad rank for min(mat)");
    assert(dim >= 0 && dim < M::order);

    Vec<Value_type<M>> result(mat.extent(dim));
    if (dim == 0) { // row
        for (std::size_t i = 0; i < mat.rows(); ++i) {
            result(i) = min(mat.row(i));
        }
    }
    else { // column
        for (std::size_t i = 0; i < mat.rows(); ++i) {
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
sum(const M& mat, std::size_t dim)
{
    static_assert(M::order == 2, "bad rank for sum(mat)");
    assert(dim >= 0 && dim < M::order);

    Vec<Value_type<M>> result(mat.extent(dim));
    if (dim == 0) { // row
        for (std::size_t i = 0; i < mat.rows(); ++i) {
            result(i) = sum(mat.row(i));
        }
    }
    else { // column
        for (std::size_t i = 0; i < mat.rows(); ++i) {
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
prod(const M& mat, std::size_t dim)
{
    static_assert(M::order == 2, "bad rank for prod(mat)");
    assert(dim >= 0 && dim < M::order);

    Vec<Value_type<M>> result(mat.extent(dim));
    if (dim == 0) { // row
        for (std::size_t i = 0; i < mat.rows(); ++i) {
            result(i) = prod(mat.row(i));
        }
    }
    else { // column
        for (std::size_t i = 0; i < mat.rows(); ++i) {
            result(i) = prod(mat.column(i));
        }
    }
    return result;
}

//------------------------------------------------------------------------------
//
// Compute trace of a square matrix:

template <typename M>
inline Enable_if<Matrix_type<M>(), typename M::value_type> trace(const M& mat)
{
    static_assert(M::order == 2, "trace: bad matrix rank");
    assert(mat.rows() == mat.cols());

    constexpr auto zero = Value_type<M>{0};

    const auto d = mat.diag();
    return std::accumulate(d.begin(), d.end(), zero);
}

//------------------------------------------------------------------------------
//
// Vector norm:
//
// TODO: Only Euclidean vector norm is implemented so far.

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

//------------------------------------------------------------------------------
//
// Vector dot and cross products:

template <typename M1, typename M2>
inline Enable_if<Matrix_type<M1>() && Matrix_type<M2>, typename M1::value_type>
dot(const M1& x, const M2& y)
{
    static_assert(M1::order == 1, "bad rank for dot product");
    static_assert(M2::order == 1, "bad rank for dot product");
    assert(same_extents(x, y));

    constexpr auto zero = Value_type<M1>{0};
    return std::inner_product(x.begin(), x.end(), y.begin(), zero);
}

template <typename T>
inline void
cross(const Matrix<T, 1>& x, const Matrix<T, 1>& y, Matrix<T, 1>& res)
{
    assert(x.size() == 3 && x.size() == y.size());
    res.resize(3);
    res(0) = x(1) * y(2) - x(2) * y(1);
    res(1) = x(2) * y(0) - x(0) * y(2);
    res(2) = x(0) * y(1) - x(1) * y(0);
}

template <typename M1, typename M2>
inline Enable_if<Matrix_type<M1>() && Matrix_type<M2>(),
                 Matrix<typename M1::value_type, 1>>
cross(const M1& x, const M2& y)
{
    Matrix<typename M1::value_type, 1> res;
    cross(x, y, res);
    return res;
}

//------------------------------------------------------------------------------

// Compute vector-scalar product and add the result to a vector.
template <typename T>
inline void axpy(const T& a, const Matrix<T, 1>& x, Matrix<T, 1>& y)
{
    assert(same_extents(x, y));
    for (std::size_t i = 0; i < x.size(); ++i) {
        y(i) = a * x(i) + y(i);
    }
}

//------------------------------------------------------------------------------

// Determinant of square matrix.
double det(const Mat<double>& a);

//------------------------------------------------------------------------------
//
// Matrix decomposition:

// LU factorization.
void lu(Mat<double>& a, Vec<int>& ipiv);

} // namespace Numlib

#endif // NUMLIB_MATH_LINALG_H
