// Copyright (c) 2018 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#include <numlib/constants.h>
#include <numlib/math.h>
#include <numlib/matrix.h>
#include <catch/catch.hpp>
#include <functional>
#include <cmath>
#include <limits>

double f(double x) { return x * x; }

TEST_CASE("test_math_calculus")
{
    SECTION("derivation")
    {
        auto dx = Numlib::Math::dfdx(f, 2.0);
        CHECK(std::abs(dx - 4.0) < 1.0e-8);
    }

    SECTION("trapz")
    {
        double xlo = 2.1;
        double xup = 3.6;

        Numlib::Vec<double> y = {3.2, 2.7, 2.9, 3.5, 4.1, 5.2};

        double ft = Numlib::Math::trapz(xlo, xup, y);

        CHECK(std::abs(ft - 5.22) < 1.0e-8);
    }

    SECTION("quad")
    {
        using namespace Numlib::Math;

        double a = 0.0;
        double b = Numlib::Constants::pi;

        double res = quad<5>([](double x) { return std::sin(x); }, a, b);
        CHECK(std::abs(res - 2.0) < 5.0e-7);

        res = quad<8>([](double x) { return std::sin(x); }, a, b);
        CHECK(std::abs(res - 2.0) < 1.0e-14);

        double eps = std::numeric_limits<double>::epsilon();
        res = quad<16>([](double x) { return std::sin(x); }, a, b);
        CHECK(std::abs(res - 2.0) < eps);
    }
}
