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
    SECTION("max_vec")
    {
        using namespace Numlib;

        Vec<int> a = {1, 2, 3, 4};
        CHECK(max(a) == 4);
    }

    SECTION("max_mat")
    {
        using namespace Numlib;

        Mat<int> m = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
        Vec<int> max_row = {3, 6, 9};
        Vec<int> max_col = {7, 8, 9};

        CHECK(max(m, 0) == max_row);
        CHECK(max(m, 1) == max_col);
    }
}
