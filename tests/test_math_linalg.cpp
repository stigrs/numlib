// Copyright (c) 2018 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#include <numlib/math.h>
#include <numlib/matrix.h>
#include <catch2/catch.hpp>

TEST_CASE("test_math_linalg")
{
    using namespace Numlib;

    SECTION("linspace")
    {
        Vec<double> ans = {0.0, 0.2, 0.4, 0.6, 0.8, 1.0};
        auto v = linspace(0.0, 1.0, 6);

        for (Index i = 0; i < v.size(); ++i) {
            CHECK(std::abs(v(i) - ans(i)) < 1.0e-12);
        }
    }

    SECTION("identity")
    {
        Mat<int> eye = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
        CHECK(identity<int>(3) == eye);
    }

    SECTION("max_min_vec")
    {
        Vec<int> a = {1, 2, 3, 4};
        CHECK(max(a) == 4);
        CHECK(min(a) == 1);
    }

    SECTION("sum_vec")
    {
        Vec<int> a = {1, 2, 3, 4};
        CHECK(sum(a) == 10);
    }

    SECTION("prod_vec")
    {
        Vec<int> a = {1, 2, 3, 4};
        CHECK(prod(a) == 24);
    }

    SECTION("max_min_mat")
    {
        Mat<int> m = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};

        Vec<int> max_row = {3, 6, 9};
        Vec<int> max_col = {7, 8, 9};

        Vec<int> min_row = {1, 4, 7};
        Vec<int> min_col = {1, 2, 3};

        CHECK(max(m, 0) == max_row);
        CHECK(max(m, 1) == max_col);
        CHECK(min(m, 0) == min_row);
        CHECK(min(m, 1) == min_col);
    }

    SECTION("sum_mat")
    {
        Mat<int> m = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};

        Vec<int> sum_row = {6, 15, 24};
        Vec<int> sum_col = {12, 15, 18};

        CHECK(sum(m, 0) == sum_row);
        CHECK(sum(m, 1) == sum_col);
    }

    SECTION("prod_mat")
    {
        Mat<int> m = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};

        Vec<int> prod_row = {6, 120, 504};
        Vec<int> prod_col = {28, 80, 162};

        CHECK(prod(m, 0) == prod_row);
        CHECK(prod(m, 1) == prod_col);
    }

    SECTION("norm_vector")
    {
        Vec<double> v = {1.0, 2.0, 3.0};
        auto vn = norm(v);
        CHECK(vn * vn == 14.0);
    }

    SECTION("norm_sparse_vector")
    {
        const double ans = 37.416573867739416; // Numpy

        Sparse_vector<double> v = {{1, 10.0}, {4, 20.0}, {9, 30.0}};
        CHECK(std::abs(norm(v) - ans) < 1.0e-12);
    }

    SECTION("norm_matrix")
    {
        // Numpy:
        Mat<double> a = {{-4.0, -3.0, -2.0}, {-1.0, 0.0, 1.0}, {2.0, 3.0, 4.0}};

        CHECK(std::abs(norm(a, 'F') - 7.745966692414834) < 1.0e-12);
        CHECK(std::abs(norm(a, 'I') - 9.0) < 1.0e-12);
        CHECK(std::abs(norm(a, '1') - 7.0) < 1.0e-12);
    }

    SECTION("normalize")
    {
        Vec<double> v = {1.0, 2.0, 3.0};
        Vec<double> vn = normalize(v);
        CHECK(vn == v / std::sqrt(14.0));
    }

    SECTION("trace")
    {
        Mat<int> a = {{-1, 0, 3}, {11, 5, 2}, {6, 12, -6}};
        const auto asub = a(slice(0, 2), slice(0, 2));

        CHECK(trace(a) == -2);
        CHECK(trace(asub) == 4);
    }

    SECTION("dot_dense_vector")
    {
        Vec<int> a = {1, 3, -5};
        Vec<int> b = {4, -2, -1};

        CHECK(dot(a, b) == 3);
    }

    SECTION("dot_sparse_vector")
    {
        const double ans = 420.0; // Numpy

        Sparse_vector<double> x = {{1, 10.0}, {4, 20.0}, {9, 30.0}};
        Vec<double> y = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0};

        CHECK(std::abs(dot(y, x) - ans) < 1.0e-12);
    }

    SECTION("cross")
    {
        Vec<double> a = {3.0, -3.0, 1.0};
        Vec<double> b = {4.0, 9.0, 2.0};
        Vec<double> axb = {-15.0, -2.0, 39.0};

        CHECK(cross(a, b) == axb);
    }

    SECTION("transpose")
    {
        Mat<int> m = {{1, 2}, {3, 4}, {5, 6}};
        Mat<int> ans = {{1, 3, 5}, {2, 4, 6}};

        CHECK(transpose(m) == ans);
    }

    SECTION("det")
    {
        const double ans2 = 13.0;
        const double ans3 = 76.0;
        const double ans4 = 242.0; // armadillo

        Mat<double> a2 = {{1.0, 5.0}, {-2.0, 3.0}};

        Mat<double> a3 = {{1.0, 5.0, 4.0}, {-2.0, 3.0, 6.0}, {5.0, 1.0, 0.0}};

        Mat<double> a4 = {{1.0, 5.0, 4.0, 2.0},
                          {-2.0, 3.0, 6.0, 4.0},
                          {5.0, 1.0, 0.0, -1.0},
                          {2.0, 3.0, -4.0, 0.0}};

        CHECK(std::abs(det(a2) - ans2) < 1.0e-12);
        CHECK(std::abs(det(a3) - ans3) < 1.0e-12);
        CHECK(std::abs(det(a4) - ans4) < 1.0e-12);
    }

    SECTION("inv")
    {
        Mat<double> a = {{1.0, 5.0, 4.0, 2.0},
                         {-2.0, 3.0, 6.0, 4.0},
                         {5.0, 1.0, 0.0, -1.0},
                         {2.0, 3.0, -4.0, 0.0}};

        // Numpy:
        Mat<double> ainv = {{-0.19008264, 0.16528926, 0.28099174, 0.05785124},
                            {0.34710744, -0.21487603, -0.16528926, 0.02479339},
                            {0.16528926, -0.0785124, 0.01652893, -0.20247934},
                            {-0.60330579, 0.61157025, 0.23966942, 0.31404959}};

        inv(a);

        for (Index i = 0; i < a.rows(); ++i) {
            for (Index j = 0; j < a.cols(); ++j) {
                CHECK(std::abs(a(i, j) - ainv(i, j)) < 1.0e-8);
            }
        }
    }

    SECTION("qr")
    {
        // Numpy:
        Mat<double> rans = {
            {-14., -21., 14.}, {0.0, -175.0, 70.0}, {0.0, 0.0, -35.0}};
        Mat<double> qans = {{-0.85714286, 0.39428571, 0.33142857},
                            {-0.42857143, -0.90285714, -0.03428571},
                            {0.28571429, -0.17142857, 0.94285714}};
        Mat<double> a = {
            {12.0, -51.0, 4.0}, {6.0, 167.0, -68.0}, {-4.0, 24.0, -41.0}};

        Mat<double> q;
        Mat<double> r;

        qr(a, q, r);

        for (Index i = 0; i < q.rows(); ++i) {
            for (Index j = 0; j < q.cols(); ++j) {
                CHECK(std::abs(q(i, j) - qans(i, j)) < 1.0e-8);
            }
        }
        for (Index i = 0; i < r.rows(); ++i) {
            for (Index j = 0; j < r.cols(); ++j) {
                CHECK(std::abs(r(i, j) - rans(i, j)) < 1.0e-12);
            }
        }

        Mat<double> qr = q * r;

        for (Index i = 0; i < a.rows(); ++i) {
            for (Index j = 0; j < a.cols(); ++j) {
                CHECK(std::abs(a(i, j) - qr(i, j)) < 1.0e-12);
            }
        }
    }

    SECTION("svd")
    {
        // Example from Intel MKL:
        Vec<double> sans = {27.47, 22.64, 8.56, 5.99, 2.01};
        Mat<double> uans = {{-0.59, 0.26, 0.36, 0.31, 0.23},
                            {-0.40, 0.24, -0.22, -0.75, -0.36},
                            {-0.03, -0.60, -0.45, 0.23, -0.31},
                            {-0.43, 0.24, -0.69, 0.33, 0.16},
                            {-0.47, -0.35, 0.39, 0.16, -0.52},
                            {0.29, 0.58, -0.02, 0.38, -0.65}};
        Mat<double> vtans = {{-0.25, -0.40, -0.69, -0.37, -0.41},
                             {0.81, 0.36, -0.25, -0.37, -0.10},
                             {-0.26, 0.70, -0.22, 0.39, -0.49},
                             {0.40, -0.45, 0.25, 0.43, -0.62},
                             {-0.22, 0.14, 0.59, -0.63, -0.44}};

        Mat<double> a = {{8.79, 9.93, 9.83, 5.45, 3.16},
                         {6.11, 6.91, 5.04, -0.27, 7.98},
                         {-9.15, -7.93, 4.86, 4.85, 3.01},
                         {9.57, 1.64, 8.83, 0.74, 5.8},
                         {-3.49, 4.02, 9.80, 10.00, 4.27},
                         {9.84, 0.15, -8.99, -6.02, -5.31}};

        Vec<double> s;
        Mat<double> u;
        Mat<double> vt;

        svd(a, s, u, vt);

        for (Index i = 0; i < sans.size(); ++i) {
            CHECK(std::abs(s(i) - sans(i)) < 5.0e-3);
        }
        for (Index i = 0; i < uans.rows(); ++i) {
            for (Index j = 0; j < uans.cols(); ++j) {
                CHECK(std::abs(u(i, j) - uans(i, j)) < 5.0e-3);
            }
        }
        for (Index i = 0; i < vtans.rows(); ++i) {
            for (Index j = 0; j < vtans.cols(); ++j) {
                CHECK(std::abs(vt(i, j) - vtans(i, j)) < 5.0e-3);
            }
        }
    }

    SECTION("eigs_dense")
    {
        // Numpy:
        Vec<double> eval = {3.28792877e-06, 3.05898040e-04, 1.14074916e-02,
                            2.08534219e-01, 1.56705069e+00};
        Mat<double> evec = {{-0.0062, 0.0472, 0.2142, -0.6019, 0.7679},
                            {0.1167, -0.4327, -0.7241, 0.2759, 0.4458},
                            {-0.5062, 0.6674, -0.1205, 0.4249, 0.3216},
                            {0.7672, 0.2330, 0.3096, 0.4439, 0.2534},
                            {-0.3762, -0.5576, 0.5652, 0.4290, 0.2098}};

        Mat<double> a = {{1.0, 0.5, 1. / 3., 1. / 4., 1. / 5},
                         {0.5, 1. / 3., 1. / 4., 1. / 5., 1. / 6.},
                         {1. / 3., 1. / 4., 1. / 5., 1. / 6., 1. / 7.},
                         {1. / 4., 1. / 5., 1. / 6., 1. / 7., 1. / 8.},
                         {1. / 5., 1. / 6., 1. / 7., 1. / 8., 1. / 9.}};

        Vec<double> w;
        eigs(a, w);

        for (int i = 0; i < 5; ++i) {
            CHECK(std::abs(w(i) - eval(i)) < 1.0e-8);
        }
        for (int i = 0; i < 5; ++i) {
            for (int j = 0; j < 5; ++j) {
                CHECK(std::abs(a(i, j) - evec(i, j)) < 1.0e-4);
            }
        }
    }

    SECTION("eig_dense")
    {
        // Numpy:
        Vec<double> eval_re = {-3.17360337, -3.17360337, 2.84219813,
                               7.50500862};

        Vec<double> eval_im = {1.12844169, -1.12844169, 0.0, 0.0};

        Mat<double> evec_re = {
            {-0.16889612, -0.16889612, -0.19514446, 0.70845976},
            {0.61501958, 0.61501958, 0.08601687, 0.46590401},
            {-0.19838031, -0.19838031, -0.58764782, 0.52110625},
            {-0.72497646, -0.72497646, 0.78050610, 0.09729590}};

        Mat<double> evec_im = {{-0.11229493, 0.11229493, 0.0, 0.0},
                               {-0.03942734, 0.03942734, 0.0, 0.0},
                               {0.11880544, -0.11880544, 0.0, 0.0},
                               {0.0, 0.0, 0.0, 0.0}};

        Mat<double> a = {{1.0, 5.0, 4.0, 2.0},
                         {-2.0, 3.0, 6.0, 4.0},
                         {5.0, 1.0, 0.0, -1.0},
                         {2.0, 3.0, -4.0, 0.0}};

        Vec<std::complex<double>> eval(4);
        Mat<std::complex<double>> evec(4, 4);

        eig(a, evec, eval);

        for (Index i = 0; i < eval.size(); ++i) {
            CHECK(std::abs(eval(i).real() - eval_re(i)) < 5.0e-8);
            CHECK(std::abs(eval(i).imag() - eval_im(i)) < 5.0e-8);
        }

        for (Index i = 0; i < evec.rows(); ++i) {
            for (Index j = 0; j < evec.cols(); ++j) {
                CHECK(std::abs(evec(i, j).real() - evec_re(i, j)) < 5.0e-9);
                CHECK(std::abs(evec(i, j).imag() - evec_im(i, j)) < 5.0e-9);
            }
        }
    }

    SECTION("eigs_band_matrix")
    {
        // Example from scipy:
        Mat<double> a = {{1.0, 5.0, 2.0, 0.0},
                         {5.0, 2.0, 5.0, 2.0},
                         {2.0, 5.0, 3.0, 5.0},
                         {0.0, 2.0, 5.0, 4.0}};

        Vec<double> w = {-4.26200532, -2.22987175, 3.95222349, 12.53965359};
        Mat<double> v = {{0.54585106, -0.73403852, 0.39896071, -0.06375287},
                         {-0.49026342, 0.04827646, 0.67105283, -0.55407514},
                         {-0.58943537, -0.41123929, 0.15802574, 0.6771086},
                         {0.33801529, 0.53827417, 0.60460427, 0.48006276}};

        Band_matrix<double> ab(3, 3, a);
        Mat<double> evec;
        Vec<double> eval;

        eigs(ab, evec, eval);

        for (Index i = 0; i < eval.size(); ++i) {
            CHECK(std::abs(eval(i) - w(i)) < 5.0e-9);
        }
        for (Index i = 0; i < evec.rows(); ++i) {
            for (Index j = 0; j < evec.cols(); ++j) {
                CHECK(std::abs(evec(i, j) - v(i, j)) < 1.0e-8);
            }
        }
    }

#ifdef USE_MKL
    SECTION("eig_band_matrix")
    {
        // Example from Intel MKL:

        double A[77] = {0.0, 0.0, 0.0, 5.0, 2.0, 1.0, 1.0, 0.0, 0.0, 2.0, 6.0,
                        3.0, 1.0, 1.0, 0.0, 1.0, 3.0, 6.0, 3.0, 1.0, 1.0, 1.0,
                        1.0, 3.0, 6.0, 3.0, 1.0, 1.0, 1.0, 1.0, 3.0, 6.0, 3.0,
                        1.0, 1.0, 1.0, 1.0, 3.0, 6.0, 3.0, 1.0, 1.0, 1.0, 1.0,
                        3.0, 6.0, 3.0, 1.0, 1.0, 1.0, 1.0, 3.0, 6.0, 3.0, 1.0,
                        1.0, 1.0, 1.0, 3.0, 6.0, 3.0, 1.0, 0.0, 1.0, 1.0, 3.0,
                        6.0, 2.0, 0.0, 0.0, 1.0, 1.0, 2.0, 5.0, 0.0, 0.0, 0.0};

        double eval_ans[6];
        eval_ans[0] = 3.1715728752538100;
        eval_ans[1] = 4.0000000000000000;
        eval_ans[2] = 4.0000000000000000;
        eval_ans[3] = 4.1292484841890931;
        eval_ans[4] = 4.4066499006731521;
        eval_ans[5] = 6.0000000000000000;

        Band_mat<double> ab(11, 11, 3, 3, A);
        Mat<double> evec(ab.rows(), ab.cols());
        Vec<double> eval(ab.cols());

        double emin = 3.0;
        double emax = 7.0;

        eig(emin, emax, ab, evec, eval);
        for (int i = 0; i < 6; ++i) {
            CHECK(std::abs(eval(i) - eval_ans[i]) < 1.0e-12);
        }
    }
#endif

#if defined(__unix__) || defined(__APPLE__) || defined(USE_MKL)
    // OpenBLAS v0.2.14.1 gives wrong results.
    SECTION("eigs_packed_matrix")
    {
        // Results from numpy:
        Vec<double> w = {3.28792877e-06, 3.05898040e-04, 1.14074916e-02,
                         2.08534219e-01, 1.56705069e+00};

        Mat<double> v = {
            {0.00617386, 0.04716181, -0.21421362, -0.60187148, -0.76785474},
            {-0.11669275, -0.43266733, 0.72410213, 0.27591342, -0.44579106},
            {0.50616366, 0.66735044, 0.12045328, 0.42487662, -0.32157829},
            {-0.76719119, 0.23302452, -0.30957397, 0.44390304, -0.25343894},
            {0.37624555, -0.55759995, -0.56519341, 0.42901335, -0.20982264}};

        auto a = hilbert<>(5);
        Symm_mat<double, lo> ap(a);

        Mat<double> evec;
        Vec<double> eval(5);

        eigs(ap, evec, eval);

        for (Index i = 0; i < eval.size(); ++i) {
            CHECK(std::abs(eval(i) - w(i)) < 1.0e-8);
        }
        for (Index i = 0; i < evec.rows(); ++i) {
            for (Index j = 0; j < evec.cols(); ++j) {
                CHECK(std::abs(evec(i, j) - v(i, j)) < 1.0e-8);
            }
        }
    }
#endif

#ifdef USE_MKL
    SECTION("eig_sparse_matrix")
    {
        // Example from Intel MKL:

        // clang-format off
        BLAS_INT rows[12] = {0, 4, 9, 15, 22, 29, 36, 43, 50, 56, 61, 65};
        BLAS_INT cols[65] = {0,   1,   2,   3,
                             0,   1,   2,   3,   4,
                             0,   1,   2,   3,   4,   5,
                             0,   1,   2,   3,   4,   5,   6,
                                  1,   2,   3,   4,   5,   6,   7,
                                       2,   3,   4,   5,   6,   7,   8,
                                            3,   4,   5,   6,   7,   8,  9,
                                                 4,   5,   6,   7,   8,  9,  10,
                                                      5,   6,   7,   8,  9,  10,
                                                           6,   7,   8,  9,  10,
                                                                7,   8,  9,  10
        };
        double val[65] = {5.0, 2.0, 1.0, 1.0,
                          2.0, 6.0, 3.0, 1.0, 1.0,
                          1.0, 3.0, 6.0, 3.0, 1.0, 1.0,
                          1.0, 1.0, 3.0, 6.0, 3.0, 1.0, 1.0,
                               1.0, 1.0, 3.0, 6.0, 3.0, 1.0, 1.0,
                                    1.0, 1.0, 3.0, 6.0, 3.0, 1.0, 1.0,
                                         1.0, 1.0, 3.0, 6.0, 3.0, 1.0, 1.0,
                                              1.0, 1.0, 3.0, 6.0, 3.0, 1.0, 1.0,
                                                   1.0, 1.0, 3.0, 6.0, 3.0, 1.0,
                                                        1.0, 1.0, 3.0, 6.0, 2.0,
                                                             1.0, 1.0, 2.0, 5.0
        };
        // clang-format on

        double eig_ans[6];
        eig_ans[0] = 3.1715728752538100;
        eig_ans[1] = 4.0000000000000000;
        eig_ans[2] = 4.0000000000000000;
        eig_ans[3] = 4.1292484841890931;
        eig_ans[4] = 4.4066499006731521;
        eig_ans[5] = 6.0000000000000000;

        Sp_mat<double> a(11, 11, val, cols, rows);
        Mat<double> evec(a.rows(), a.cols());
        Vec<double> eval(a.cols());

        double emin = 3.0;
        double emax = 7.0;

        eig(emin, emax, a, evec, eval);

        for (int i = 0; i < 6; ++i) {
            CHECK(std::abs(eval(i) - eig_ans[i]) < 1.0e-12);
        }
    }
#endif

    SECTION("linsolve")
    {
        Mat<double> A = {{1.0, 2.0, 3.0}, {2.0, 3.0, 4.0}, {3.0, 4.0, 1.0}};
        Mat<double> B = {{14.0}, {20.0}, {14.0}};
        Vec<double> x = {1.0, 2.0, 3.0};

        linsolve(A, B);

        for (Index i = 0; i < B.rows(); ++i) {
            CHECK(std::abs(B(i, 0) - x(i)) < 1.0e-12);
        }
    }

#ifdef USE_MKL
    SECTION("sparse_linsolve")
    {
        // Example from Matlab:
        Vec<double> xans = {1.0, 2.0, 3.0, 4.0, 5.0};

        Mat<double> A = {{0.0, 2.0, 0.0, 1.0, 0.0},
                         {4.0, -1.0, -1.0, 0.0, 0.0},
                         {0.0, 0.0, 0.0, 3.0, -6.0},
                         {-2.0, 0.0, 0.0, 0.0, 2.0},
                         {0.0, 0.0, 4.0, 2.0, 0.0}};

        Sp_mat<double> SA = gather(A);

        Mat<double> B = {{8.0}, {-1.0}, {-18.0}, {8.0}, {20.0}};
        Mat<double> x;

        linsolve(SA, B, x);

        for (Index i = 0; i < xans.size(); ++i) {
            CHECK(std::abs(x(i, 0) - xans(i)) < 1.0e-12);
        }
    }
#endif

    SECTION("rcond")
    {
        Mat<double> A = hilbert(10);
        Mat<double> B = identity(5) * 0.01;

        CHECK(std::abs(rcond(A) - 2.8286e-14) < 1.0e-12);
        CHECK(std::abs(rcond(B) - 1.0) < 1.0e-12);
    }

    SECTION("kron")
    {
        Mat<int> ans = {{1, 1, 0, 0}, {1, 1, 0, 0}, {0, 0, 1, 1}, {0, 0, 1, 1}};
        Mat<int> res;
        kron(identity<int>(2), ones<Mat<int>>(2, 2), res);

        CHECK(res == ans);
    }

    SECTION("lstsq")
    {
        // Example from Intel MKL:
        Mat<double> xans = {{-0.69, -0.24, 0.06},
                            {-0.80, -0.08, 0.21},
                            {0.38, 0.12, -0.65},
                            {0.29, -0.24, 0.42},
                            {0.29, 0.35, -0.30}};

        Mat<double> a = {{0.12, -8.19, 7.69, -2.26, -4.71},
                         {-6.91, 2.22, -5.12, -9.08, 9.96},
                         {-3.33, -8.94, -6.72, -4.40, -9.98},
                         {3.97, 3.33, -2.74, -7.92, -3.20}};

        Mat<double> b = {{7.30, 0.47, -6.28},
                         {1.33, 6.58, -3.42},
                         {2.68, -1.71, 3.46},
                         {-9.62, -0.79, 0.41},
                         {0.00, 0.00, 0.00}};

        lstsq(a, b);

        for (Index i = 0; i < b.rows(); ++i) {
            for (Index j = 0; j < b.cols(); ++j) {
                CHECK(std::abs(b(i, j) - xans(i, j)) < 5.0e-3);
            }
        }

        // Example from Numpy:

        Mat<double> an = {{0.0, 1.0}, {1.0, 1.0}, {2.0, 1.0}, {3.0, 1.0}};
        Mat<double> y = {{-1.0}, {0.2}, {0.9}, {2.1}};

        lstsq(an, y);

        CHECK(std::abs(y(0, 0) - 1.0) < 1.0e-12);
        CHECK(std::abs(y(1, 0) + 0.95) < 1.0e-12);
    }
}
