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

    SECTION("norm")
    {
        Vec<double> v = {1.0, 2.0, 3.0};
        auto vn = norm(v);
        CHECK(vn * vn == 14.0);
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

    SECTION("dot")
    {
        Vec<int> a = {1, 3, -5};
        Vec<int> b = {4, -2, -1};

        CHECK(dot(a, b) == 3);
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

    SECTION("eigs")
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

    SECTION("eig")
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
}
