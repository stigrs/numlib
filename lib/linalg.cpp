// Copyright (c) 2018 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#include <numlib/math.h>
#include <cmath>

Numlib::Vec<double> Numlib::linspace(double x1, double x2, Index n)
{
    Vec<double> res(n);
    Vec<double>::iterator it;

    double h = (x2 - x1) / (n - 1);
    double val = 0.0;

    for (it = res.begin(), val = x1; it != res.end(); ++it, val += h) {
        *it = val;
    }
    return res;
}

double Numlib::det(const Mat<double>& a)
{
    assert(a.rows() == a.cols());

    double ddet = 0.0;
    const Index n = a.rows();

    if (n == 1) {
        ddet = a(0, 0);
    }
    else if (n == 2) {
        ddet = a(0, 0) * a(1, 1) - a(1, 0) * a(0, 1);
    }
    else { // use LU decomposition
        Mat<double> tmp(a);
        Vec<Index> ipiv;

        lu(tmp, ipiv);

        Index permut = 0;
        for (Index i = 1; i <= n; ++i) {
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

    const Index n = a.cols();

    evec.resize(n, n);
    eval.resize(n);

    Vec<double> wr(n);
    Vec<double> wi(n);
    Mat<double> vr(n, n);
    Mat<double> vl(n, n);

    Index info =
        LAPACKE_dgeev(LAPACK_ROW_MAJOR, 'N', 'V', n, a.data(), n, wr.data(),
                      wi.data(), vl.data(), n, vr.data(), n);
    if (info != 0) {
        throw Math_error("dgeev failed");
    }
    for (Index i = 0; i < vr.rows(); ++i) {
        std::complex<double> wii(wr(i), wi(i));
        eval(i) = wii;
        for (Index j = 0; j < vr.cols(); j += 2) {
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

#ifdef USE_MKL
void Numlib::eig(double emin,
                 double emax,
                 const Numlib::Band_mat<double>& ab,
                 Numlib::Mat<double>& evec,
                 Numlib::Vec<double>& eval)
{
    // Intitialize FEAST:

    Index fpm[128];
    feastinit((Index*) fpm);
#ifndef NDEBUG
    fpm[0] = 1; // print runtime status
#endif

    // Solve eigenvalue problem:

    Index n = ab.cols();
    Index kla = ab.upper();
    Index lda = ab.leading_dim();
    Index m0 = n;
    Index loop = 0;
    Index m = m0;
    Index info = 0;

    double epsout = 0.0; // relative error on the trace (not returned)
    Vec<double> res(m0); // residual vector (not returned)

    evec.resize(n, m0);
    eval.resize(m0);

    dfeast_sbev("F", &n, &kla, ab.data(), &lda, (Index*) fpm, &epsout, &loop,
                &emin, &emax, &m0, eval.data(), evec.data(), &m, res.data(),
                &info);
    if (info != 0) {
        throw Math_error("dfeast_sbev failed");
    }

    // Return the m first eigenvalues and eigenvectors:

    eval = eval(slice{0, m});
    evec = evec(slice{0, n}, slice{0, m});
}
#endif

#ifdef USE_MKL
void Numlib::eig(double emin,
                 double emax,
                 const Numlib::Sp_mat<double>& a,
                 Numlib::Mat<double>& evec,
                 Numlib::Vec<double>& eval)
{
    // Initialize FEAST:

    Index fpm[128];
    feastinit((Index*) fpm);
#ifndef NDEBUG
    fpm[0] = 1; // print runtime status
#endif

    // Solve eigenvalue problem:

    Index n = a.cols();
    Index m0 = n;
    Index loop = 0;
    Index m = m0;
    Index info = 0;

    auto ia = a.row_index_one_based(); // FEAST only support one-based indexing
    auto ja = a.columns_one_based();

    double epsout = 0.0; // relative error on the trace (not returned)
    Vec<double> res(m0); // residual vector (not returned)

    evec.resize(n, m0);
    eval.resize(m0);

    dfeast_scsrev("F", &n, a.data(), ia.data(), ja.data(), (Index*) fpm,
                  &epsout, &loop, &emin, &emax, &m0, eval.data(), evec.data(),
                  &m, res.data(), &info);
    if (info != 0) {
        throw Math_error("dfeast_scsrev failed");
    }

    // Return the m first eigenvalues and eigenvectors:

    eval = eval(slice{0, m});
    evec = evec(slice{0, n}, slice{0, m});
}
#endif

void Numlib::schmidt(Mat<double>& a, Index n)
{
    Index n_out = 0;
    Index n_orb = n;
    Index n_bas = a.rows();

    Vec<double> work(n_bas);
    work = 0.0;

    double r_min = 0.1;

    while (n_orb < n_bas) {
        Index lim = n_orb + n_bas;
        for (Index i = 0; i < lim; ++i) {
            if (n_out >= n_bas) {
                return;
            }
            auto an = a.column(n_out);
            if (i < n_orb) {
                auto ai = a.column(i);
                an = ai;
            }
            else {
                an = 0.0;
                a(i - n_orb, n_out) = 1.0;
            }
            for (Index j = 0; j < n_out; ++j) {
                auto aj = a.column(j);
                work(j) = dot(aj, an);
            }
            for (Index j = 0; j < n_out; ++j) {
                auto aj = a.column(j);
                an = an - work(j) * aj;
            }
            double r = std::sqrt(dot(an, an));
            if (r >= r_min) {
                ++n_out;
                an /= r;
            }
        }
        r_min /= 10.0;
        n_orb = n_out;
    }
}
