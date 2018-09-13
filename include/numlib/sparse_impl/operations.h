// Copyright (c) 2018 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#ifndef NUMLIB_SPARSE_OPERATIONS_H
#define NUMLIB_SPARSE_OPERATIONS_H

#include <numlib/matrix.h>
#include <vector>
#include <iostream>

namespace Numlib {

//------------------------------------------------------------------------------
//
// Convert formats:

// Gather a sparse full-storage vector into compressed form.
template <typename T>
Sparse_vector<T> gather(const Matrix<T, 1>& y)
{
    std::vector<T> val;
    std::vector<std::ptrdiff_t> loc;

    for (std::ptrdiff_t i = 0; i < y.size(); ++i) {
        if (y(i) != T{0}) {
            val.push_back(y(i));
            loc.push_back(i);
        }
    }
    return {val, loc};
}

// Gather a sparse full-storage matrix into sparse CSR3 format.
template <typename T>
Sparse_matrix<T> gather(const Matrix<T, 2>& m)
{
    std::vector<T> values;
    std::vector<std::ptrdiff_t> columns;
    std::vector<std::ptrdiff_t> row_index(m.rows() + 1);

    std::ptrdiff_t nnz = 0;
    for (std::ptrdiff_t i = 0; i < m.rows(); ++i) {
        std::ptrdiff_t inz = 0;
        for (std::ptrdiff_t j = 0; j < m.cols(); ++j) {
            if (m(i, j) != T{0}) {
                values.push_back(m(i, j));
                columns.push_back(j);
                nnz++;
                inz++;
            }
        }
        row_index[i] = nnz - inz;
    }
    row_index[row_index.size() - 1] = nnz;
    return {m.rows(), m.cols(), values, columns, row_index};
}

// Scatter a sparse vector into full storage form.
template <typename T>
Matrix<T, 1> scatter(const Sparse_vector<T>& y)
{
    Matrix<T, 1> res(y.size());
    for (std::ptrdiff_t i = 0; i < res.size(); ++i) {
        res(i) = y(i);
    }
    return res;
}

// Scatter a sparse matrix into full storage form.
template <typename T>
Matrix<T, 2> scatter(const Sparse_matrix<T>& m)
{
    Matrix<T, 2> res(m.rows(), m.cols());

    for (std::ptrdiff_t i = 0; i < res.rows(); ++i) {
        for (std::ptrdiff_t j = 0; j < res.cols(); ++j) {
            res(i, j) = m(i, j);
        }
    }
    return res;
}

//------------------------------------------------------------------------------
//
// Binary arithmetic operations:

// Scalar multiplication:

template <typename T>
inline Sparse_vector<T> operator*(const Sparse_vector<T>& a, const T& scalar)
{
    Sparse_vector<T> res(a);
    return res *= scalar;
}

template <typename T>
inline Sparse_vector<T> operator*(const T& scalar, const Sparse_vector<T>& a)
{
    Sparse_vector<T> res(a);
    return res *= scalar;
}

template <typename T>
inline Sparse_matrix<T> operator*(const Sparse_matrix<T>& a, const T& scalar)
{
    Sparse_matrix<T> res(a);
    return res *= scalar;
}

template <typename T>
inline Sparse_matrix<T> operator*(const T& scalar, const Sparse_matrix<T>& a)
{
    Sparse_matrix<T> res(a);
    return res *= scalar;
}

// Scalar division:

template <typename T>
inline Sparse_vector<T> operator/(const Sparse_vector<T>& a, const T& scalar)
{
    Sparse_vector<T> res(a);
    return res /= scalar;
}

template <typename T>
inline Sparse_matrix<T> operator/(const Sparse_matrix<T>& a, const T& scalar)
{
    Sparse_matrix<T> res(a);
    return res /= scalar;
}

// Vector addition:

template <typename T>
Matrix<T, 1> operator+(const Sparse_vector<T>& x, const Matrix<T, 1>& y)
{
    assert(x.size() == y.size());

    Matrix<T, 1> res(y);

    std::ptrdiff_t i = 0;
    for (const auto& v : x) {
        res(x.loc(i)) += v;
        ++i;
    }
    return res;
}

template <typename T>
Matrix<T, 1> operator+(const Matrix<T, 1>& y, const Sparse_vector<T>& x)
{
    assert(x.size() == y.size());

    Matrix<T, 1> res(y);

    std::ptrdiff_t i = 0;
    for (const auto& v : x) {
        res(x.loc(i)) += v;
        ++i;
    }
    return res;
}

// Vector subtraction:

template <typename T>
Matrix<T, 1> operator-(const Sparse_vector<T>& x, const Matrix<T, 1>& y)
{
    assert(x.size() == y.size());

    Matrix<T, 1> res(y);

    std::ptrdiff_t i = 0;
    for (const auto& v : x) {
        res(x.loc(i)) -= v;
        ++i;
    }
    return res;
}

template <typename T>
Matrix<T, 1> operator-(const Matrix<T, 1>& y, const Sparse_vector<T>& x)
{
    assert(x.size() == y.size());

    Matrix<T, 1> res(y);

    std::ptrdiff_t i = 0;
    for (const auto& v : x) {
        res(x.loc(i)) -= v;
        ++i;
    }
    return res;
}

//------------------------------------------------------------------------------
//
// Matrix-vector product:

template <typename T>
void mv_mul(const Sparse_matrix<T>& a, const Matrix<T, 1>& x, Matrix<T, 1>& res)
{
    assert(x.size() == a.cols());

    using size_type = typename Matrix<T, 1>::size_type;

    res.resize(a.cols());

    for (size_type i = 0; i < a.rows(); ++i) {
        T sum = T{0};
        for (size_type j = a.row_index()[i]; j < a.row_index()[i + 1]; ++j) {
            sum += a.values()[j] * x(a.columns()[j]);
        }
        res(i) = sum;
    }
}

template <typename T>
inline Matrix<T, 1> operator*(const Sparse_matrix<T>& a, const Matrix<T, 1>& x)
{
    assert(x.size() == a.cols());

    Matrix<T, 1> res(a.cols());
    mv_mul(a, x, res);
    return res;
}

//------------------------------------------------------------------------------
//
// Output to stream:

// Output stream operator for sparse vectors.
template <typename T>
std::ostream& operator<<(std::ostream& to, const Sparse_vector<T>& vec)
{
    std::ptrdiff_t i = 0;

    to << "[number of non-zero elements: " << vec.num_nonzero() << "]\n";
    for (const auto& x : vec) {
        to << "(" << vec.loc(i) << ")\t" << x << '\n';
        ++i;
    }
    to << '\n';
    return to;
}

// Output stream operator for sparse matrices.
template <class T>
std::ostream& operator<<(std::ostream& to, const Sparse_matrix<T>& mat)
{
    to << "[matrix size: " << mat.rows() << " x " << mat.cols()
       << "; number of non-zero elements: " << mat.num_nonzero() << "]\n\n";
    for (std::ptrdiff_t i = 0; i < mat.rows(); ++i) {
        for (std::ptrdiff_t j = 0; j < mat.cols(); ++j) {
            if (mat(i, j) != T{0}) {
                to << "(" << i << ", " << j << ")\t" << mat(i, j) << '\n';
            }
        }
    }
    to << '\n';
    return to;
}

} // namespace Numlib

#endif // NUMLIB_SPARSE_OPERATIONS_H
