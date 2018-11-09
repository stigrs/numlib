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

    // Answers are computed with WolframAlpha:

    SECTION("comp_ellint_1")
    {
        CHECK(std::abs(comp_ellint_1(0.0) - 1.570796326794896619) < 1.0e-15);
        CHECK(std::abs(comp_ellint_1(0.5) - 1.685750354812596043) < 1.0e-15);
        CHECK(std::abs(comp_ellint_1(0.9) - 2.280549138422770205) < 1.0e-15);
    }

    SECTION("comp_ellint_2")
    {
        CHECK(std::abs(comp_ellint_2(0.0) - 1.570796326794896619) < 1.0e-15);
        CHECK(std::abs(comp_ellint_2(0.5) - 1.467462209339427155) < 1.0e-15);
        CHECK(std::abs(comp_ellint_2(0.9) - 1.171697052781614141) < 1.0e-15);
    }
}

