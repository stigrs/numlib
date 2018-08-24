// Copyright (c) 2018 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#include <numlib/math.h>
#include <numlib/matrix.h>
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

    SECTION("integration")
    {
        double xlo = 2.1;
        double xup = 3.6;

        Numlib::Vec<double> y = {3.2, 2.7, 2.9, 3.5, 4.1, 5.2};

        double ft = Numlib::Math::trapezoidal(xlo, xup, y);

        CHECK(std::abs(ft - 5.22) < 1.0e-8);
    }
}
