// Copyright (c) 2018 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#include <numlib/matrix.h>
#include <catch/catch.hpp>
#include <iostream>

TEST_CASE("test_matrix2")
{
    using namespace Numlib;

    Matrix<int, 2> m2 = {{1, 2, 3, 4}, {5, 6, 7, 8}, {9, 10, 11, 12}};

    SECTION("rank") { CHECK(m2.rank() == 2); }
    SECTION("size") { CHECK(m2.size() == 12); }

    SECTION("extents")
    {
        CHECK(m2.extent(0) == 3);
        CHECK(m2.extent(1) == 4);
    }

    SECTION("subscripting")
    {
        int it = 1;
        for (std::size_t i = 0; i < m2.extent(0); ++i) {
            for (std::size_t j = 0; j < m2.extent(1); ++j) {
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

        for (std::size_t i = 0; i < r0.size(); ++i) {
            CHECK(r0(i) == m2_r0(i));
        }
        for (std::size_t i = 0; i < r1.size(); ++i) {
            CHECK(r1(i) == m2_r1(i));
        }
        for (std::size_t i = 0; i < r2.size(); ++i) {
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

        for (std::size_t i = 0; i < c0.size(); ++i) {
            CHECK(c0(i) == m2_c0(i));
        }
        for (std::size_t i = 0; i < c1.size(); ++i) {
            CHECK(c1(i) == m2_c1(i));
        }
        for (std::size_t i = 0; i < c2.size(); ++i) {
            CHECK(c2(i) == m2_c2(i));
        }
        for (std::size_t i = 0; i < c3.size(); ++i) {
            CHECK(c3(i) == m2_c3(i));
        }
    }

    SECTION("slice")
    {
        Matrix<int, 2> m3 = {{01, 02, 03}, {11, 12, 13}, {21, 22, 23}};

        auto m30 = m3(Slice{0, 2}, Slice{0, 2});
        auto m31 = m3(Slice{1, 2}, 1);
        auto m32 = m3(Slice{1, 2}, 0);

        for (auto x : m30) {
            std::cout << x << std::endl;
        }
        CHECK(m30.rank() == 2);
        CHECK(m30(0, 0) == 01);
        CHECK(m30(1, 0) == 11);

        CHECK(m31.rank() == 2);
        CHECK(m31(0, 0) == 12);
        CHECK(m31(1, 0) == 22);

        CHECK(m32.rank() == 2);
        CHECK(m32(0, 0) == 11);
        CHECK(m32(1, 0) == 21);
    }
}
