// Copyright (c) 2018 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#include <numlib/math.h>
#include <catch2/catch.hpp>
#include <cmath>

TEST_CASE("test_math_special_functions")
{
    using namespace Numlib;

    SECTION("comp_ellint_1")
    {
        CHECK(std::abs(comp_ellint_1(0.0) - 1.5707963267948966) < 5.12e-15);
        CHECK(std::abs(comp_ellint_1(0.5) - 1.6857503548125961) < 5.12e-15);
        CHECK(std::abs(comp_ellint_1(0.9) - 2.2805491384227703) < 5.12e-15);
    }

    SECTION("comp_ellint_2")
    {
        CHECK(std::abs(comp_ellint_2(0.0) - 1.5707963267948966) < 2.0e-9);
        CHECK(std::abs(comp_ellint_2(0.5) - 1.4674622093394272) < 2.0e-9);
        CHECK(std::abs(comp_ellint_2(0.9) - 1.1716970527816142) < 2.0e-9);
    }
}

