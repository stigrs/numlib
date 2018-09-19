// Copyright (c) 2018 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#include <numlib/math.h>
#include <numlib/matrix.h>
#include <catch2/catch.hpp>

TEST_CASE("test_math_euler")
{
    using namespace Numlib;

    SECTION("eul2rotm_identity")
    {
        Mat<double> rotm = eul2rotm();
        Mat<double> eye = identity(3);

        CHECK(same_extents(rotm, eye));

        for (Index i = 0; i < rotm.rows(); ++i) {
            for (Index j = 0; j < rotm.cols(); ++j) {
                CHECK(std::abs(rotm(i, j) - eye(i, j)) < 1.0e-12);
            }
        }
    }

    SECTION("eul2rotm")
    {
        Mat<double> ans = {{0.0, 0.0, 1.0}, {0.0, 1.0, 0.0}, {-1.0, 0.0, 0.0}};
        Mat<double> rotm = eul2rotm(0.0, 90.0, 0.0);

        CHECK(same_extents(rotm, ans));

        for (Index i = 0; i < rotm.rows(); ++i) {
            for (Index j = 0; j < rotm.cols(); ++j) {
                CHECK(std::abs(rotm(i, j) - ans(i, j)) < 1.0e-12);
            }
        }
    }

    SECTION("rotm2eul_identity")
    {
        Vec<double> ans = {0.0, 0.0, 0.0};
        Vec<double> res = rotm2eul(identity(3));

        CHECK(same_extents(res, ans));

        for (Index i = 0; i < res.size(); ++i) {
            CHECK(std::abs(res(i) - ans(i)) < 1.0e-12);
        }
    }

    SECTION("rotm2eul")
    {
        Mat<double> rotm = {{0.0, 0.0, 1.0}, {0.0, 1.0, 0.0}, {-1.0, 0.0, 0.0}};
        Vec<double> ans = {0.0, 90.0, 0.0};
        Vec<double> res = rotm2eul(rotm);

        CHECK(same_extents(res, ans));

        for (Index i = 0; i < res.size(); ++i) {
            CHECK(std::abs(res(i) - ans(i)) < 1.0e-12);
        }
    }

    SECTION("eul2quat")
    {
        Vec<double> ans = {0.7071, 0.0, 0.7071, 0.0};
        Vec<double> res = eul2quat(0.0, 90.0, 0.0);

        CHECK(same_extents(res, ans));

        for (Index i = 0; i < res.size(); ++i) {
            CHECK(std::abs(res(i) - ans(i)) < 1.0e-5);
        }
    }
}
