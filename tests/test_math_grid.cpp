// Copyright (c) 2018 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#include <numlib/math.h>
#include <catch2/catch.hpp>

TEST_CASE("test_math_geometry")
{
    using namespace Numlib;

    SECTION("2d-meshgrid")
    {
        Mat<double> ans = {{2.0, 5.0, 10.0},
                           {5.0, 8.0, 13.0},
                           {10.0, 13.0, 18.0},
                           {17.0, 20.0, 25.0},
                           {26.0, 29.0, 34.0}};

        Vec<double> x = {1.0, 2.0, 3.0};
        Vec<double> y = {1.0, 2.0, 3.0, 4.0, 5.0};

        Mat<double> xx;
        Mat<double> yy;

        meshgrid(x, y, xx, yy);

        auto zz = pow(xx, 2.0) + pow(yy, 2.0);

        for (Index i = 0; i < zz.rows(); ++i) {
            for (Index j = 0; j < zz.cols(); ++j) {
                CHECK(zz(i, j) == ans(i, j));
            }
        }
    }

    SECTION("3d-meshgrid")
    {
        Mat<double> ans0 = {{0, 4, 16, 36},  {1, 5, 17, 37},   {4, 8, 20, 40},
                            {9, 13, 25, 45}, {16, 20, 32, 52}, {25, 29, 41, 61},
                            {36, 40, 52, 72}};
        Mat<double> ans1 = {{9, 13, 25, 45},  {10, 14, 26, 46},
                            {13, 17, 29, 49}, {18, 22, 34, 54},
                            {25, 29, 41, 61}, {34, 38, 50, 70},
                            {45, 49, 61, 81}};
        Mat<double> ans2 = {{36, 40, 52, 72}, {37, 41, 53, 73},
                            {40, 44, 56, 76}, {45, 49, 61, 81},
                            {52, 56, 68, 88}, {61, 65, 77, 97},
                            {72, 76, 88, 108}};

        Vec<double> x = {0.0, 2.0, 4.0, 6.0};
        Vec<double> y = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
        Vec<double> z = {0.0, 3.0, 6.0};

        Cube<double> xx;
        Cube<double> yy;
        Cube<double> zz;

        meshgrid(x, y, z, xx, yy, zz);

        auto f = pow(xx, 2.0) + pow(yy, 2.0) + pow(zz, 2.0);

        for (Index i = 0; i < f.extent(0); ++i) {
            for (Index j = 0; j < f.extent(1); ++j) {
                CHECK(f(i, j, 0) == ans0(i, j));
            }
        }
        for (Index i = 0; i < f.extent(0); ++i) {
            for (Index j = 0; j < f.extent(1); ++j) {
                CHECK(f(i, j, 1) == ans1(i, j));
            }
        }
        for (Index i = 0; i < f.extent(0); ++i) {
            for (Index j = 0; j < f.extent(1); ++j) {
                CHECK(f(i, j, 2) == ans2(i, j));
            }
        }
    }
}
