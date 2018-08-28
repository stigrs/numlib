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

//------------------------------------------------------------------------------
//
// Vector dot and cross products:

template <typename M>
inline Enable_if<Matrix_type<M>(), typename M::value_type> dot(const M& x,
                                                               const M& y)
{
    static_assert(M::order == 1, "bad rank for dot product");
    assert(same_extents(x, y));

    constexpr auto zero = Value_type<M>{0};
    return std::inner_product(x.begin(), x.end(), y.begin(), zero);
}

} // namespace Numlib

#endif // NUMLIB_MATH_LINALG_H
