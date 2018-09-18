// Copyright (c) 2018 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#include <numlib/matrix.h>
#include <catch2/catch.hpp>

TEST_CASE("test_matrix2")
{
    using namespace Numlib;

    Matrix<int, 2> m2 = {{1, 2, 3, 4}, {5, 6, 7, 8}, {9, 10, 11, 12}};

    SECTION("order") { CHECK(m2.rank() == 2); }
    SECTION("size") { CHECK(m2.size() == 12); }

    SECTION("extents")
    {
        CHECK(m2.extent(0) == 3);
        CHECK(m2.extent(1) == 4);
        CHECK(rows(m2) == 3);
        CHECK(cols(m2) == 4);
    }

    SECTION("subscripting")
    {
        int it = 1;
        for (Index i = 0; i < m2.extent(0); ++i) {
            for (Index j = 0; j < m2.extent(1); ++j) {
                CHECK(m2(i, j) == it);
                ++it;
            }
        }
    }

    SECTION("row")
    {
        Matrix<int, 1> m2_r0 = {1, 2, 3, 4};
        Matrix<int, 1> m2_r1 = {5, 6, 7, 8};
        Matrix<int, 1> m2_r2 = {9, 10, 11, 12};

        auto r0 = m2.row(0);
        auto r1 = m2.row(1);
        auto r2 = m2.row(2);

        for (Index i = 0; i < r0.size(); ++i) {
            CHECK(r0(i) == m2_r0(i));
        }
        for (Index i = 0; i < r1.size(); ++i) {
            CHECK(r1(i) == m2_r1(i));
        }
        for (Index i = 0; i < r2.size(); ++i) {
            CHECK(r2(i) == m2_r2(i));
        }
    }

    SECTION("column")
    {
        Matrix<int, 1> m2_c0 = {1, 5, 9};
        Matrix<int, 1> m2_c1 = {2, 6, 10};
        Matrix<int, 1> m2_c2 = {3, 7, 11};
        Matrix<int, 1> m2_c3 = {4, 8, 12};

        auto c0 = m2.column(0);
        auto c1 = m2.column(1);
        auto c2 = m2.column(2);
        auto c3 = m2.column(3);

        for (Index i = 0; i < c0.size(); ++i) {
            CHECK(c0(i) == m2_c0(i));
        }
        for (Index i = 0; i < c1.size(); ++i) {
            CHECK(c1(i) == m2_c1(i));
        }
        for (Index i = 0; i < c2.size(); ++i) {
            CHECK(c2(i) == m2_c2(i));
        }
        for (Index i = 0; i < c3.size(); ++i) {
            CHECK(c3(i) == m2_c3(i));
        }
    }

    SECTION("slice")
    {
        Matrix<int, 2> m3 = {{01, 02, 03}, {11, 12, 13}, {21, 22, 23}};

        auto m30 = m3(slice{0, 2}, slice{0, 2});
        auto m31 = m3(slice{1, 2}, 1);
        auto m32 = m3(slice{1, 2}, 0);

        CHECK(m30.rank() == 2);
        CHECK(m30.size() == 4);
        CHECK(m30(0, 0) == 01);
        CHECK(m30(1, 0) == 11);
        CHECK(m30(0, 1) == 02);
        CHECK(m30(1, 1) == 12);

        CHECK(m31.rank() == 2);
        CHECK(m31(0, 0) == 12);
        CHECK(m31(1, 0) == 22);

        CHECK(m32.rank() == 2);
        CHECK(m32(0, 0) == 11);
        CHECK(m32(1, 0) == 21);

        Matrix<int, 2> m4_ans = {{00, 00, 03}, {00, 00, 13}, {21, 22, 23}};

        Matrix<int, 2> m4(m3);
        auto m40 = m4(slice{0, 2}, slice{0, 2});

        m40 = 0;

        CHECK(m4 == m4_ans);
    }

    SECTION("construct_from_array")
    {
        int ptr[4] = {10, 20, 30, 40};

        Matrix<int, 2> m5 = Matrix_ref<int, 2>(Matrix_slice<2>{2, 2}, &ptr[0]);

        CHECK(m5(0, 0) == 10);
        CHECK(m5(0, 1) == 20);
        CHECK(m5(1, 0) == 30);
        CHECK(m5(1, 1) == 40);
    }

    SECTION("diag")
    {
        Matrix<int, 1> m6_ans = {1, 1, 1};

        auto m6 = zeros<Matrix<int, 2>>(3, 3);
        auto d = m6.diag();

        d = 1;

        CHECK(d == m6_ans);
    }

    SECTION("slice_diag")
    {
        Matrix<int, 1> s7_ans = {1, 6};

        // clang-format off
        Matrix<int, 2> m7 = {{ 1,  2,  3,  4},
                             { 5,  6,  7,  8},
                             { 9, 10, 11, 12},
                             {13, 14, 15, 16}};
        // clang-format on

        auto s7 = m7(slice{0, 2}, slice{0, 2});
        CHECK(s7.diag() == s7_ans);
    }

    SECTION("a_plus_b")
    {
        Matrix<int, 2> a = {{1, 2}, {3, 4}};
        Matrix<int, 2> b = {{10, 20}, {30, 40}};

        auto c = a + b;

        CHECK(c(0, 0) == 11);
        CHECK(c(0, 1) == 22);
        CHECK(c(1, 0) == 33);
        CHECK(c(1, 1) == 44);
    }

    SECTION("a_minus_b")
    {
        Matrix<int, 2> a = {{10, 20}, {30, 40}};
        Matrix<int, 2> b = {{1, 2}, {3, 4}};

        auto c = a - b;

        CHECK(c(0, 0) == 9);
        CHECK(c(0, 1) == 18);
        CHECK(c(1, 0) == 27);
        CHECK(c(1, 1) == 36);
    }

    SECTION("copy_ctor")
    {
        Matrix<double, 2> a = {{1.0, 2.0}, {3.0, 4.0}};
        auto b(a);
        CHECK(a == b);
        CHECK(a.size() == b.size());
        CHECK(a.rows() == b.rows());
        CHECK(a.cols() == b.cols());
    }

    SECTION("assignment")
    {
        Matrix<double, 2> a = {{1.0, 2.0}, {3.0, 4.0}};
        auto b = a;
        CHECK(a == b);
        CHECK(a.size() == b.size());
        CHECK(a.rows() == b.rows());
        CHECK(a.cols() == b.cols());
    }

    SECTION("swap")
    {
        std::swap(m2, m2);
        CHECK(m2 == m2);

        Matrix<int, 2> a1 = {{-1, 0, -6}, {6, 5, 2}, {11, 12, 3}};
        Matrix<int, 2> a2 = {{-1, 0, -6}, {6, 5, 2}, {11, 12, 3}};
        Matrix<int, 2> a3 = {{11, 12, 3}, {6, 5, 2}, {-1, 0, -6}};
        std::swap(a2, a3);
        CHECK(a3 == a1);
    }

    SECTION("swap_rows")
    {
        Matrix<int, 2> a1 = {{-1, 0, -6}, {6, 5, 2}, {11, 12, 3}};
        Matrix<int, 2> a2 = {{6, 5, 2}, {-1, 0, -6}, {11, 12, 3}};

        a1.swap_rows(0, 1);

        CHECK(a1 == a2);
    }

    SECTION("mm_mul")
    {
        Matrix<int, 2> a = {{1, 2, 3}, {4, 5, 6}};
        Matrix<int, 2> b = {{7, 8}, {9, 10}, {11, 12}};
        Matrix<int, 2> ans = {{58, 64}, {139, 154}};

        CHECK((a * b) == ans);
    }

    SECTION("mv_mul")
    {
        Matrix<int, 2> a = {{1, -1, 2}, {0, -3, 1}};
        Matrix<int, 1> x = {2, 1, 0};
        Matrix<int, 1> y = {1, -3};

        CHECK((a * x) == y);
    }
}
