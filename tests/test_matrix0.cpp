// Copyright (c) 2018 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#include <numlib/matrix.h>
#include <catch2/catch.hpp>

TEST_CASE("test_matrix0")
{
    using namespace Numlib;

    Matrix<int, 0> m0 = {1};

    SECTION("order") { CHECK(m0.order == 0); }
    SECTION("subscripting") { CHECK(m0() == 1); }
}
