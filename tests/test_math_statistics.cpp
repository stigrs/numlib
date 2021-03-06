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

    SECTION("movmean")
    {
        Vec<double> ans1 = {4.0,  8.0,  6.0, -1.0, -2.0,
                            -3.0, -1.0, 3.0, 4.0,  5.0};
        Vec<double> ans3 = {6.0,  6.0,     4.3333, 1.0, -2.0,
                            -2.0, -0.3333, 2.0,    4.0, 4.5};
        Vec<double> ans10 = {3.0, 2.0,    1.5714, 1.75,   2.0,
                             2.3, 2.1111, 1.375,  0.7143, 1.0};
        Vec<double> ans40 = {4.0, 6.0,  6.0,  4.25, 3.0,
                             1.6, -0.2, -0.8, 0.2,  1.6};

        Vec<double> v = {4.0, 8.0, 6.0, -1.0, -2.0, -3.0, -1.0, 3.0, 4.0, 5.0};

        auto res = movmean(v, 1);
        for (Index i = 0; i < v.size(); ++i) {
            CHECK(std::abs(res(i) - ans1(i)) < 5.0e-5);
        }
        res = movmean(v, 3);
        for (Index i = 0; i < v.size(); ++i) {
            CHECK(std::abs(res(i) - ans3(i)) < 5.0e-5);
        }
        res = movmean(v, 10);
        for (Index i = 0; i < v.size(); ++i) {
            CHECK(std::abs(res(i) - ans10(i)) < 5.0e-5);
        }
        res = movmean(v, {4, 0});
        for (Index i = 0; i < v.size(); ++i) {
            CHECK(std::abs(res(i) - ans40(i)) < 5.0e-5);
        }
    }

    SECTION("cov")
    {
        Vec<double> b = {3.0,  13.0, 7.0,  5.0,  21.0, 23.0, 23.0,
                         40.0, 23.0, 14.0, 12.0, 56.0, 23.0, 29.0};
        Vec<double> c = {3.0,  13.0, 7.0,  5.0,  21.0, 23.0, 39.0,
                         23.0, 40.0, 23.0, 14.0, 12.0, 56.0, 23.0};

        CHECK(std::abs(cov(b, c) - 59.78021978) < 1.0e-8);
    }

    SECTION("kabsch_rmsd")
    {
        // Example from Matlab:
        Mat<double> p = {{0.8147, 0.9058, 0.1270}, {0.9134, 0.6324, 0.0975},
                         {0.2785, 0.5469, 0.9575}, {0.9649, 0.1576, 0.9706},
                         {0.9572, 0.4854, 0.8003}, {0.1419, 0.4218, 0.9157}};

        Mat<double> q = {{0.7922, 0.9595, 0.6557}, {0.0357, 0.8491, 0.9340},
                         {0.6787, 0.7577, 0.7431}, {0.3922, 0.6555, 0.1712},
                         {0.7060, 0.0318, 0.2769}, {0.0462, 0.0971, 0.8235}};

        CHECK(std::abs(kabsch_rmsd(p, q) - 0.4761) < 5.0e-6);
    }
}
