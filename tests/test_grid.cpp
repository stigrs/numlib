// Copyright (c) 2018 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#include <numlib/math.h>
#include <catch2/catch.hpp>
#include <iostream>

TEST_CASE("test_grid")
{
    Numlib::Grid g(-1.0, 10.0, 0.1);
    CHECK(g.start() == -1.0);
    CHECK(g.max() == 10.0);
    CHECK(g.step() == 0.1);
    CHECK(g.size() == 111);
}

