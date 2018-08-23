// Copyright (c) 2018 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#include <numlib/constants.h>
#include <catch/catch.hpp>

TEST_CASE("test_constants")
{
    using namespace num;

    SECTION("constants")
    {
        CHECK(pi == 3.14159265358979323846);
        CHECK(h_bar == 1.05457180000e-34);
    }
}
