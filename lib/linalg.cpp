// Copyright (c) 2018 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#include <numlib/matrix.h>
#include <numlib/math.h>
#include <cmath>

#ifdef __APPLE__
#include <Accelerate/Accelerate.h>
#else
#include <lapacke.h>
#endif

double Numlib::det(const Mat<double>& a)
{
    assert(a.rows() == a.cols());

    double ddet = 0.0;
    const int n = static_cast<int>(a.rows());

    if (n == 1) {
        ddet = a(0, 0);
    }
    else if (n == 2) {
        ddet = a(0, 0) * a(1, 1) - a(1, 0) * a(0, 1);
    }
    else { // use LU decomposition
        Mat<double> tmp(a);
        Vec<int> ipiv;

        lu(tmp, ipiv);

        int permut = 0;
        for (int i = 1; i <= n; ++i) {
            if (i != ipiv(i - 1)) { // Fortran uses base 1
                permut++;
            }
        }
        ddet = prod(tmp.diag());
        ddet *= std::pow(-1.0, static_cast<double>(permut));
    }
    return ddet;
}

void Numlib::inv(Mat<double>& a)
{
    assert(a.rows() == a.cols());

    if (det(a) == 0.0) {
        throw Math_error("inv: matrix not invertible");
    }
    const int n = static_cast<int>(a.rows());
    const int lda = n;

    Vec<int> ipiv(n);
    lu(a, ipiv); // perform LU factorization

    int info = LAPACKE_dgetri(LAPACK_ROW_MAJOR, n, a.data(), lda, ipiv.data());
    if (info != 0) {
        throw Math_error("dgetri: matrix inversion failed");
    }
}

void Numlib::lu(Mat<double>& a, Vec<int>& ipiv)
{
    const int m = static_cast<int>(a.rows());
    const int n = static_cast<int>(a.cols());
    const int lda = n;

    ipiv.resize(std::min(m, n));

    int info =
        LAPACKE_dgetrf(LAPACK_ROW_MAJOR, m, n, a.data(), lda, ipiv.data());
    if (info < 0) {
        throw Math_error("dgetrf: illegal input parameter");
    }
    if (info > 0) {
        throw Math_error("dgetrf: U matrix is singular");
    }
}
