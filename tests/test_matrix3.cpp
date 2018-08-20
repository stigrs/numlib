// Copyright (c) 2018 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#include <numlib/matrix.h>
#include <catch/catch.hpp>
#include <iostream>

TEST_CASE("test_matrix3")
{
    using namespace Numlib;

    Matrix<int, 3> m3 = {
        {{1, 2}, {3, 4}}, {{5, 6}, {7, 8}}, {{9, 10}, {11, 12}}};

    SECTION("rank") { CHECK(m3.rank() == 3); }
    SECTION("size") { CHECK(m3.size() == 12); }

    SECTION("extents")
    {
        CHECK(m3.extent(0) == 3);
        CHECK(m3.extent(1) == 2);
        CHECK(m3.extent(2) == 2);
    }

    SECTION("subscripting")
    {
        CHECK(m3(0, 0, 0) == 1);
        CHECK(m3(0, 0, 1) == 2);
        CHECK(m3(0, 1, 0) == 3);
        CHECK(m3(0, 1, 1) == 4);
        CHECK(m3(1, 0, 0) == 5);
        CHECK(m3(1, 0, 1) == 6);
        CHECK(m3(1, 1, 0) == 7);
        CHECK(m3(1, 1, 1) == 8);
    }

    SECTION("row")
    {
        auto r = m3.row(1);
        auto rr = r.row(1);
        int rc = rr.row(1);
        CHECK(r.rank() == 2);
        CHECK(rr.rank() == 1);
        CHECK(rc == 8);
    }
}
