// Copyright (c) 2018 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#include <numlib/constants.h>
#include <numlib/math.h>
#include <numlib/matrix.h>
#include <catch2/catch.hpp>

double moda(int k, int l, int m)
{
    return narrow_cast<double>(narrow_cast<long long>(k * l) % m);
}

TEST_CASE("test_math_signal")
{
    using namespace Numlib;

    SECTION("convolution")
    {
        Vec<int> ans = {4, 13, 28, 34, 32, 21}; // numpy
        Vec<int> a = {1, 2, 3};
        Vec<int> b = {4, 5, 6, 7};
        Vec<int> c = conv(a, b);
        CHECK(c == ans);
    }

#ifdef USE_MKL
    SECTION("fft")
    {
        // Example from Intel MKL (basic_dp_real_dft_1d.c):

        const double twopi = 2.0 * Constants::pi;
        const int n = 6;
        int h = -1;

        Vec<double> x(n);

        double factor = (2 * (n - h) % n == 0) ? 1.0 : 2.0;
        for (int i = 0; i < n; ++i) {
            double phase = moda(i, h, n) / n;
            x(i) = factor * std::cos(twopi * phase) / n;
        }

        auto y = fft(x);

        const double errthr = 2.5 * std::log(narrow_cast<double>(n)) /
                              std::log(2.0) * DBL_EPSILON;

        for (int i = 0; i < y.size(); ++i) {
            double re_exp = 0.0;
            if ((i - h) % n == 0 || (-i - h) % n == 0) {
                re_exp = 1.0;
            }
            double im_exp = 0.0;
            double re_got = y(i).real();
            double im_got = y(i).imag();
            double err = std::abs(re_got - re_exp) + std::abs(im_got - im_exp);
            CHECK(err <= errthr);
        }

        auto xback = ifft(y);
        for (int i = 0; i < n; ++i) {
            double err = std::abs(xback(i) - x(i));
            CHECK(err <= errthr);
        }
    }
#endif
}
