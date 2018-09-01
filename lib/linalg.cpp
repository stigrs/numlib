// Copyright (c) 2018 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#include <numlib/math.h>
#include <cmath>

double Numlib::det(const Mat<double>& a)
{
    assert(a.rows() == a.cols());

    double ddet = 0.0;
    const int n = narrow_cast<int>(a.rows());

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
        ddet *= std::pow(-1.0, narrow_cast<double>(permut));
    }
    return ddet;
}

void Numlib::eig(Mat<double>& a,
                 Mat<std::complex<double>>& evec,
                 Vec<std::complex<double>>& eval)
{
    assert(a.rows() == a.cols());

    const int n = narrow_cast<int>(a.cols());

    evec.resize(n, n);
    eval.resize(n);

    Vec<double> wr(n);
    Vec<double> wi(n);
    Mat<double> vr(n, n);
    Mat<double> vl(n, n);

    int info = LAPACKE_dgeev(LAPACK_ROW_MAJOR, 'N', 'V', n, a.data(), n,
                             wr.data(), wi.data(), vl.data(), n, vr.data(), n);
    if (info != 0) {
        throw Math_error("dgeev failed");
    }
    for (std::size_t i = 0; i < vr.rows(); ++i) {
        std::complex<double> wii(wr(i), wi(i));
        eval(i) = wii;
        for (std::size_t j = 0; j < vr.cols(); j += 2) {
            std::complex<double> v1 = {vr(i, j), 0.0};
            std::complex<double> v2 = {vr(i, j + 1), 0.0};
            if (wi(j) != 0.0) {
                v1 = {vr(i, j), vr(i, j + 1)};
                v2 = {vr(i, j), -vr(i, j + 1)};
            }
            evec(i, j) = v1;
            evec(i, j + 1) = v2;
        }
    }
}
