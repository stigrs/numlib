// Copyright (c) 2018 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#include <numlib/math.h>
#include <catch/catch.hpp>
#include <functional>

double f(double x) { return x * x; }

TEST_CASE("test_math_calculus")
{
    SECTION("derivation")
    {
        auto dx = Numlib::Math::dfdx(f, 2.0);
        CHECK(std::abs(dx - 4.0) < 1.0e-8);
    }
}
