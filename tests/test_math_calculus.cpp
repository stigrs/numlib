// Copyright (c) 2018 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#include <numlib/constants.h>
#include <numlib/math.h>
#include <numlib/matrix.h>
#include <catch2/catch.hpp>
#include <functional>
#include <cmath>
#include <limits>

double f(double x) { return x * x; }
double sinf(double& x) { return std::sin(x); }
double invexp(double& x) { return std::exp(-x); }
double rate(double t, double y) { return t * std::sqrt(y); }

void lorenz(double /* t */,
            const Numlib::Vec<double>& y,
            Numlib::Vec<double>& ydot)
{
    const double sigma = 10.0;
    const double R = 28.0;
    const double b = 8.0 / 3.0;

    ydot(0) = sigma * (y(1) - y(0));
    ydot(1) = R * y(0) - y(1) - y(0) * y(2);
    ydot(2) = -b * y(2) + y(0) * y(1);
}

TEST_CASE("test_math_calculus")
{
    SECTION("derivation")
    {
        auto dx = Numlib::dfdx(f, 2.0);
        CHECK(std::abs(dx - 4.0) < 1.0e-8);
    }

    SECTION("trapz")
    {
        double xlo = 2.1;
        double xup = 3.6;

        Numlib::Vec<double> y = {3.2, 2.7, 2.9, 3.5, 4.1, 5.2};

        double ft = Numlib::trapz(xlo, xup, y);

        CHECK(std::abs(ft - 5.22) < 1.0e-8);
    }

    SECTION("quad")
    {
        using namespace Numlib;

        double a = 0.0;
        double b = Constants::pi;

        double res = quad<5>([](double x) { return std::sin(x); }, a, b);
        CHECK(std::abs(res - 2.0) < 5.0e-7);

        res = quad<8>([](double x) { return std::sin(x); }, a, b);
        CHECK(std::abs(res - 2.0) < 1.0e-14);

        double eps = std::numeric_limits<double>::epsilon();
        res = quad<16>([](double x) { return std::sin(x); }, a, b);
        CHECK(std::abs(res - 2.0) < eps);
    }

#ifdef ENABLE_QUADPACK
    SECTION("qags")
    {
        using namespace Numlib;

        double a = 0.0;
        double b = Constants::pi;

        double res = qags(sinf, a, b);
        double eps = std::numeric_limits<double>::epsilon();
        CHECK(std::abs(res - 2.0) < eps);
    }

    SECTION("qagi")
    {
        using namespace Numlib;

        double bound = 0.0;
        int inf = 1; // +infinity

        double res = qagi(invexp, bound, inf);
        CHECK(std::abs(res - 1.0) < 1.0e-15);
    }
#endif // ENABLE_QUADPACK

    SECTION("rk4")
    {
        double t0 = 0.0;
        double t1 = 10.0;
        double dt = 0.1;
        int n = 1 + narrow_cast<int>((t1 - t0) / dt);

        Numlib::Vec<double> y(n);
        y(0) = 1.0;

        for (int i = 1; i < n; ++i) {
            y(i) = Numlib::rk4(rate, y(i - 1), t0 + dt * (i - 1), dt);
        }

        for (int i = 0; i < n; ++i) {
            double t = t0 + dt * i;
            double y2 = std::pow(t * t + 4.0, 2.0) / 16.0;
            double err = std::abs(y(i) / y2 - 1.0);
            CHECK(err <= 3.0e-7);
        }
    }

    SECTION("dopri5")
    {
        using namespace Numlib;

        Mat<double> ans = {
            // results from Matlab
            {12.420121076782189, 22.132678932307815, 11.996473826705991},
            {19.500081683089384, 16.224736836476261, 45.258556702999961},
            {6.613599319856808, -7.931580903108999, 37.735650643710017},
            {-2.963989264539828, -8.250556890143775, 28.287476810924446},
            {-6.217033890199554, -8.278471219613175, 25.168552598624345}};

        Vec<double> y = {10.0, 1.0, 1.0};

        double t0 = 0.0;
        double t1 = 0.1;

        for (int i = 0; i < 5; ++i) {
            dopri5(lorenz, y, t0, t1);
            for (int j = 0; j < y.size(); ++j) {
                CHECK(std::abs(y(j) - ans(i, j)) < 1.0e-5);
            }
            t1 += 0.1;
        }
    }
}
