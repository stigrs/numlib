// Copyright (c) 2018 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#include <numlib/matrix.h>
#include <catch/catch.hpp>
#include <iostream>

TEST_CASE("test_matrix1")
{
    using namespace num;

    ivec m1 = {1, 2, 3, 4};

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
        Matrix<int, 1> m2(m1);
        auto s = m2(slice{0, 3});

        s = 0;

        CHECK(s.rank() == 1);
        CHECK(s.size() == 3);
        CHECK(s(0) == 0);
        CHECK(s(1) == 0);
        CHECK(s(2) == 0);
        CHECK(m2(0) == 0);
        CHECK(m2(1) == 0);
        CHECK(m2(2) == 0);
        CHECK(m2(3) == 4);
        CHECK(m1(0) == 1);
        CHECK(m1(1) == 2);
        CHECK(m1(2) == 3);
        CHECK(m1(3) == 4);

        auto ss = s(slice{1, 2});
        ss = 1;

        CHECK(s(0) == 0);
        CHECK(s(1) == 1);
        CHECK(s(2) == 1);
        CHECK(m2(0) == 0);
        CHECK(m2(1) == 1);
        CHECK(m2(2) == 1);
        CHECK(m2(3) == 4);
        CHECK(m1(0) == 1);
        CHECK(m1(1) == 2);
        CHECK(m1(2) == 3);
        CHECK(m1(3) == 4);
    }

    SECTION("head")
    {
        auto h = m1(slice{0, 3});

        CHECK(h.size() == 3);
        CHECK(h(0) == 1);
        CHECK(h(1) == 2);
        CHECK(h(2) == 3);
    }

    SECTION("tail")
    {
        auto t = m1(slice{2});

        CHECK(t.size() == 2);
        CHECK(t(0) == 3);
        CHECK(t(1) == 4);
    }

    SECTION("slice_stride_2")
    {
        auto s = m1(slice{0, 2, 2});

        CHECK(s.size() == 2);
        CHECK(s(0) == 1);
        CHECK(s(1) == 3);
    }

    SECTION("vector_addition")
    {
        ivec v1 = {2, 4, 5};
        ivec v2 = {2, 4, 6};
        ivec ans = {4, 8, 11};

        CHECK((v1 + v2) == ans);
    }
}
