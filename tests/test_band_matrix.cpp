// Copyright (c) 2018 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#include <numlib/matrix.h>
#include <numlib/band_matrix.h>
#include <catch/catch.hpp>

TEST_CASE("test_band_matrix")
{
    SECTION("element_access")
    {
        using namespace Numlib;

        // clang-format off
        Mat<int> a = {{11, 12,  0,  0,  0},
                      {21, 22, 23,  0,  0},
                      {31, 32, 33, 34,  0},
                      { 0, 42, 43, 44, 45},
                      { 0,  0, 53, 54, 55}};
        // clang-format on

        Band_matrix<int> ab(2, 1, a);
#
        CHECK(ab(0, 0) == 11);
        CHECK(ab(0, 1) == 12);
        CHECK(ab(1, 0) == 21);
        CHECK(ab(1, 1) == 22);
        CHECK(ab(1, 2) == 23);
        CHECK(ab(2, 0) == 31);
        CHECK(ab(2, 1) == 32);
        CHECK(ab(2, 2) == 33);
        CHECK(ab(2, 3) == 34);
        CHECK(ab(3, 1) == 42);
        CHECK(ab(3, 2) == 43);
        CHECK(ab(3, 3) == 44);
        CHECK(ab(3, 4) == 45);
        CHECK(ab(4, 2) == 53);
        CHECK(ab(4, 3) == 54);
        CHECK(ab(4, 4) == 55);
    }
}
