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
        CHECK(std::abs(comp_ellint_1(0.0) - 1.5707963267948966) < 2.0e-8);
        CHECK(std::abs(comp_ellint_1(0.5) - 1.8540746773013719) < 2.0e-8);
        CHECK(std::abs(comp_ellint_1(0.9) - 2.5780921133481733) < 2.0e-8);
    }

    SECTION("comp_ellint_2")
    {
        CHECK(std::abs(comp_ellint_2(0.0) - 1.5707963267948966) < 2.0e-8);
        CHECK(std::abs(comp_ellint_2(0.5) - 1.3506438810476755) < 2.0e-8);
        CHECK(std::abs(comp_ellint_2(0.9) - 1.1047747327040733) < 2.0e-8);
    }
}

