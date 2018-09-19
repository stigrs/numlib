// Copyright (c) 2018 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#include <numlib/math.h>
#include <numlib/matrix.h>
#include <catch2/catch.hpp>

TEST_CASE("test_math_statistics")
{
    using namespace Numlib;

    SECTION("mean_median_stddev_rms")
    {
        Vec<double> a = {3.0,  13.0, 7.0,  5.0,  21.0, 23.0, 39.0, 23.0,
                         40.0, 23.0, 14.0, 12.0, 56.0, 23.0, 29.0};

        CHECK(std::abs(mean(a) - 22.066666666666666) < 1.0e-8);
        CHECK(std::abs(median(a) - 23.0) < 1.0e-8);
        CHECK(std::abs(stddev(a) - 14.49860420211283) < 1.0e-8);
        CHECK(std::abs(rms(a) - 26.136819495365792) < 1.0e-8);
    }

    SECTION("cov")
    {
        Vec<double> b = {3.0,  13.0, 7.0,  5.0,  21.0, 23.0, 23.0,
                         40.0, 23.0, 14.0, 12.0, 56.0, 23.0, 29.0};
        Vec<double> c = {3.0,  13.0, 7.0,  5.0,  21.0, 23.0, 39.0,
                         23.0, 40.0, 23.0, 14.0, 12.0, 56.0, 23.0};

        CHECK(std::abs(cov(b, c) - 59.78021978) < 1.0e-8);
    }
}
