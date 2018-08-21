// Copyright (c) 2018 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#include <numlib/matrix.h>
#include <catch/catch.hpp>
#include <iostream>

TEST_CASE("test_matrix4")
{
    using namespace num;

    auto m4 = zeros<Matrix<int, 4>>(2, 3, 4, 5);

    SECTION("rank") { CHECK(m4.rank() == 4); }
    SECTION("size") { CHECK(m4.size() == 120); }

    SECTION("extents")
    {
        CHECK(m4.extent(0) == 2);
        CHECK(m4.extent(1) == 3);
        CHECK(m4.extent(2) == 4);
        CHECK(m4.extent(3) == 5);
    }

    SECTION("subscripting")
    {
        m4(0, 0, 0, 0) = 1;
        m4(0, 0, 0, 1) = 2;

        CHECK(m4(0, 0, 0, 0) == 1);
        CHECK(m4(0, 0, 0, 1) == 2);
    }
}
