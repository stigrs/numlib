// Copyright (c) 2018 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#include <numlib/math.h>
#include <numlib/matrix.h>
#include <catch2/catch.hpp>

TEST_CASE("test_math_quaternion")
{
    using namespace Numlib;

    SECTION("quat2rotm")
    {
        Vec<double> quat = {0.707107, 0.0, 0.707107, 0.0};
        Mat<double> ans = {{0.0, 0.0, 1.0}, {0.0, 1.0, 0.0}, {-1.0, 0.0, 0.0}};
        Mat<double> rotm = quat2rotm(quat);

        CHECK(same_extents(rotm, ans));

        for (Index i = 0; i < rotm.rows(); ++i) {
            for (Index j = 0; j < rotm.cols(); ++j) {
                CHECK(std::abs(rotm(i, j) - ans(i, j)) < 1.0e-12);
            }
        }
    }
}
