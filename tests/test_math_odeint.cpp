// Copyright (c) 2018 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#include <numlib/math.h>
#include <catch2/catch.hpp>

void fsys(double /* t */, double* y, double* ydot, void* /* data */)
{
    ydot[0] = -0.04 * y[0] + 1.0e4 * y[1] * y[2];
    ydot[2] = 3.0e7 * y[1] * y[1];
    ydot[1] = -ydot[0] - ydot[2];
}

TEST_CASE("test_math_odeint")
{
    using namespace Numlib;

    SECTION("robertson")
    {
        Mat<double> ans = {{9.851712e-01, 3.386380e-05, 1.479493e-02},
                           {9.055333e-01, 2.240655e-05, 9.444430e-02},
                           {7.158403e-01, 9.186334e-06, 2.841505e-01},
                           {4.505250e-01, 3.222964e-06, 5.494717e-01},
                           {1.831976e-01, 8.941773e-07, 8.168015e-01},
                           {3.898729e-02, 1.621940e-07, 9.610125e-01},
                           {4.936362e-03, 1.984221e-08, 9.950636e-01},
                           {5.161833e-04, 2.065787e-09, 9.994838e-01},
                           {5.179804e-05, 2.072027e-10, 9.999482e-01},
                           {5.283675e-06, 2.113481e-11, 9.999947e-01},
                           {4.658667e-07, 1.863468e-12, 9.999995e-01},
                           {1.431100e-08, 5.724404e-14, 1.000000e+00}};

        Vec<double> rtol = {1.0e-4, 1.0e-8, 1.0e-4};
        Vec<double> atol = {1.0e-6, 1.0e-10, 1.0e-6};

        Odeint ode(fsys, rtol, atol);

        Vec<double> y = {1.0, 0.0, 0.0};

        double t0 = 0.0;
        double t1 = 0.4;

        for (int i = 0; i < 12; ++i) {
            ode.integrate(y, t0, t1);
            for (Index j = 0; j < y.size(); ++j) {
                CHECK(std::abs(y(j) - ans(i, j)) / ans(i, j) < 5.0e-4);
            }
            t1 *= 10.0;
        }
    }
}
