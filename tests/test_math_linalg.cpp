// Copyright (c) 2018 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#include <numlib/math.h>
#include <numlib/matrix.h>
#include <catch/catch.hpp>

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
}
