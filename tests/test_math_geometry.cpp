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

    SECTION("pdist_matrix")
    {
        Mat<double> mat(4, 3);
        mat(0, 0) = 0.0;
        mat(0, 1) = 0.0;
        mat(0, 2) = 0.0;
        mat(1, 0) = 1.0;
        mat(1, 1) = 1.0;
        mat(1, 2) = 1.0;
        mat(2, 0) = 2.0;
        mat(2, 1) = 2.0;
        mat(2, 2) = 2.0;
        mat(3, 0) = 3.0;
        mat(3, 1) = 3.0;
        mat(3, 2) = 3.0;

        Mat<double> dm;
        pdist_matrix(dm, mat);

        Mat<double> dm_ans(4, 4);
        dm_ans = 0.0;

        dm_ans(0, 1) = 1.73205081;
        dm_ans(0, 2) = 3.46410162;
        dm_ans(0, 3) = 5.19615242;
        dm_ans(1, 2) = 1.73205081;
        dm_ans(1, 3) = 3.46410162;
        dm_ans(2, 3) = 1.73205081;
        dm_ans(1, 0) = dm_ans(0, 1);
        dm_ans(2, 0) = dm_ans(0, 2);
        dm_ans(3, 0) = dm_ans(0, 3);
        dm_ans(2, 1) = dm_ans(1, 2);
        dm_ans(3, 1) = dm_ans(1, 3);
        dm_ans(3, 2) = dm_ans(2, 3);

        for (Index i = 0; i < dm.rows(); ++i) {
            for (Index j = 0; j < dm.cols(); ++j) {
                CHECK(std::abs(dm(i, j) - dm_ans(i, j)) < 5.0e-9);
            }
        }
    }

    SECTION("meshgrid")
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
}
