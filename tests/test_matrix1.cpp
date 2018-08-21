// Copyright (c) 2018 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#include <numlib/matrix.h>
#include <catch/catch.hpp>

TEST_CASE("test_matrix1")
{
    using namespace num;

    Matrix<int, 1> m1 = {1, 2, 3, 4};

    SECTION("rank") { CHECK(m1.rank() == 1); }
    SECTION("size") { CHECK(m1.size() == 4); }

    SECTION("subscripting")
    {
        CHECK(m1(0) == 1);
        CHECK(m1(1) == 2);
        CHECK(m1(2) == 3);
        CHECK(m1(3) == 4);
    }

    SECTION("row")
    {
        int r = m1.row(1);
        CHECK(r == 2);
    }

    SECTION("slice")
    {
        auto s = m1(Slice{0, 3});

        CHECK(s.rank() == 1);
        CHECK(s.size() == 3);
        CHECK(s(0) == 1);
        CHECK(s(1) == 2);
        CHECK(s(2) == 3);
    }
}
