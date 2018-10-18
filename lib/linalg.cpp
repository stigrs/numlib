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

void Numlib::schmidt(Mat<double>& a, Index n)
{
    Index n_out = 0;
    Index n_orb = n;
    Index n_bas = a.cols();

    Vec<double> work(n_bas);
    work = 0.0;

    double r_min = 0.1;

    while (n_orb < n_bas) {
        Index lim = n_orb + n_bas;
        for (Index i = 0; i < lim; ++i) {
            if (n_out >= n_bas) {
                return;
            }
            auto an = a.row(n_out);
            if (i < n_orb) {
                auto ai = a.row(i);
                an = ai;
            }
            else {
                an = 0.0;
                a(n_out, i - n_orb) = 1.0;
            }
            for (Index j = 0; j < n_out; ++j) {
                auto aj = a.row(j);
                work(j) = dot(aj, an);
            }
            for (Index j = 0; j < n_out; ++j) {
                auto aj = a.row(j);
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

void Numlib::eig(Mat<double>& a,
                 Mat<std::complex<double>>& evec,
                 Vec<std::complex<double>>& eval)
{
    assert(a.rows() == a.cols());

    const BLAS_INT n = narrow_cast<BLAS_INT>(a.cols());

    evec.resize(n, n);
    eval.resize(n);

    Vec<double> wr(n);
    Vec<double> wi(n);
    Mat<double> vr(n, n);
    Mat<double> vl(n, n);

    BLAS_INT info =
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

    MKL_INT fpm[128];
    feastinit((MKL_INT*) fpm);
#ifndef NDEBUG
    fpm[0] = 1; // print runtime status
#endif

    // Solve eigenvalue problem:

    MKL_INT n = narrow_cast<MKL_INT>(ab.cols());
    MKL_INT kla = narrow_cast<MKL_INT>(ab.upper());
    MKL_INT lda = narrow_cast<MKL_INT>(ab.leading_dim());
    MKL_INT m0 = n;
    MKL_INT loop = 0;
    MKL_INT m = m0;
    MKL_INT info = 0;

    double epsout = 0.0; // relative error on the trace (not returned)
    Vec<double> res(m0); // residual vector (not returned)

    evec.resize(n, m0);
    eval.resize(m0);

    dfeast_sbev("F", &n, &kla, ab.data(), &lda, (MKL_INT*) fpm, &epsout, &loop,
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

    MKL_INT fpm[128];
    feastinit((MKL_INT*) fpm);
#ifndef NDEBUG
    fpm[0] = 1; // print runtime status
#endif

    // Solve eigenvalue problem:

    MKL_INT n = narrow_cast<MKL_INT>(a.cols());
    MKL_INT m0 = n;
    MKL_INT loop = 0;
    MKL_INT m = m0;
    MKL_INT info = 0;

    auto ia = a.row_index_one_based(); // FEAST only support one-based indexing
    auto ja = a.columns_one_based();

    double epsout = 0.0; // relative error on the trace (not returned)
    Vec<double> res(m0); // residual vector (not returned)

    evec.resize(n, m0);
    eval.resize(m0);

    dfeast_scsrev("F", &n, a.data(), ia.data(), ja.data(), (MKL_INT*) fpm,
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

double Numlib::det(const Mat<double>& a)
{
    assert(a.rows() == a.cols());

    double ddet = 0.0;
    const BLAS_INT n = narrow_cast<BLAS_INT>(a.rows());

    if (n == 1) {
        ddet = a(0, 0);
    }
    else if (n == 2) {
        ddet = a(0, 0) * a(1, 1) - a(1, 0) * a(0, 1);
    }
    else { // use LU decomposition
        Mat<double> tmp(a);
        Vec<BLAS_INT> ipiv;

        lu(tmp, ipiv);

        BLAS_INT permut = 0;
        for (BLAS_INT i = 1; i <= n; ++i) {
            if (i != ipiv(i - 1)) { // Fortran uses base 1
                permut++;
            }
        }
        ddet = prod(tmp.diag());
        ddet *= std::pow(-1.0, narrow_cast<double>(permut));
    }
    return ddet;
}

#ifdef USE_MKL
void Numlib::linsolve(const Numlib::Sp_mat<double>& a,
                      Numlib::Mat<double>& b,
                      Numlib::Mat<double>& x)
{
    assert(b.rows() == a.rows());
    x.resize(b.rows(), b.cols());

    MKL_INT n = narrow_cast<MKL_INT>(b.rows());
    MKL_INT nrhs = narrow_cast<MKL_INT>(b.cols());

    // Initialize PARDISO:

    void* pt[64];       // internal solver memory pointer
    MKL_INT iparm[64];  // PARDISO control parameters
    MKL_INT mtype = 11; // real and nonsymmetric matrix
    MKL_INT maxfct = 1; // max factors kept in memory
    MKL_INT mnum = 1;   // which matrix to factorize
    MKL_INT phase = 13; // analysis, numerical factorization, solve
    MKL_INT msglvl = 0; // no print of statistical information
    MKL_INT error = 0;  // initialize error flag

    pardisoinit((void*) pt, &mtype, (MKL_INT*) iparm); // set default values
    iparm[34] = 1;                                     // zero-based indexing

    // Solve linear system of equations:

    Numlib::Vec<MKL_INT> perm(n);

    pardiso((void*) pt, &maxfct, &mnum, &mtype, &phase, &n, a.data(),
            a.row_index().data(), a.columns().data(), perm.data(), &nrhs,
            (MKL_INT*) iparm, &msglvl, b.data(), x.data(), &error);

    if (error != 0) {
        throw Math_error("could not solve sparse linear system of equations");
    }
}
#endif
