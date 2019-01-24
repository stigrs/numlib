// Copyright (c) 2018 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#include <numlib/math.h>
#include <numlib/constants.h>
#include <catch2/catch.hpp>

TEST_CASE("test_math_core")
{
    using namespace Numlib;

    SECTION("even_odd")
    {
        int even_number = 4;
        CHECK(even(even_number) == true);

        int odd_number = -3;
        CHECK(odd(odd_number) == true);

        unsigned even_unumber = 4;
        CHECK(even(even_unumber) == true);

        unsigned odd_unumber = 3;
        CHECK(odd(odd_unumber) == true);
    }

    SECTION("krond")
    {
        int i = 2;
        int j = 3;
        CHECK(krond(i, i) == 1);
        CHECK(krond(i, j) == 0);

        unsigned ui = 2;
        unsigned uj = 3;
        CHECK(krond(ui, ui) == 1);
        CHECK(krond(ui, uj) == 0);
    }

    SECTION("round")
    {
        double x = 2.4;
        CHECK(round<int>(x) == 2);

        x = 2.9;
        CHECK(round<int>(x) == 3);
    }

    SECTION("square")
    {
        int x = 2;
        CHECK(square(x) == 4);
    }

    SECTION("cube")
    {
        int x = 2;
        CHECK(cube(x) == 8);
    }

    SECTION("abs")
    {
        Vec<int> v = {-1, -2, -3};
        Vec<int> ans = {1, 2, 3};
        auto res = abs(v);
        CHECK(res == ans);
    }

    SECTION("pow")
    {
        Vec<double> v = {1.0, 2.0, 3.0};
        Vec<double> ans = {1.0, 4.0, 9.0};
        auto res = pow(v, 2.0);
        for (Index i = 0; i < v.size(); ++i) {
            CHECK(res(i) == ans(i));
        }
    }

    SECTION("sqrt")
    {
        Vec<double> ans = {1.0, 2.0, 3.0};
        Vec<double> v = {1.0, 4.0, 9.0};
        auto res = sqrt(v);
        for (Index i = 0; i < v.size(); ++i) {
            CHECK(res(i) == ans(i));
        }
    }

    SECTION("sin")
    {
        using namespace Constants;

        Vec<double> ans = {0.0, 0.0, 0.0};
        Vec<double> v = {pi, pi, pi};
        auto res = sin(v);
        for (Index i = 0; i < v.size(); ++i) {
            CHECK(std::abs(res(i) - ans(i)) < 1.0e-15);
        }
    }

    SECTION("cos")
    {
        using namespace Constants;

        Vec<double> ans = {-1.0, -1.0, -1.0};
        Vec<double> v = {pi, pi, pi};
        auto res = cos(v);
        for (Index i = 0; i < v.size(); ++i) {
            CHECK(std::abs(res(i) - ans(i)) < 1.0e-15);
        }
    }

    SECTION("conj")
    {
        Mat<std::complex<double>> ans = {{{1.0, -0.0}, {1.0, -1.0}},
                                         {{-2.0, 1.0}, {0.0, -1.0}},
                                         {{5.0, -0.0}, {4.0, 2.0}}};

        Mat<std::complex<double>> m = {{{1.0, 0.0}, {1.0, 1.0}},
                                       {{-2.0, -1.0}, {0.0, 1.0}},
                                       {{5.0, 0.0}, {4.0, -2.0}}};

        auto res = conj(m);
        for (Index i = 0; i < res.rows(); ++i) {
            for (Index j = 0; j < res.cols(); ++j) {
                CHECK(res(i, j) == ans(i, j));
            }
        }
    }
}
