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
    using namespace Numlib;

    Matrix<int, 1> m1 = {1, 2, 3, 4};

    SECTION("size") { CHECK(m1.size() == 4); }

    SECTION("subscripting")
    {
        CHECK(m1(0) == 1);
        CHECK(m1(1) == 2);
        CHECK(m1(2) == 3);
        CHECK(m1(3) == 4);
    }
}
