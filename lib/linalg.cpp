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
