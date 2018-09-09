// Copyright (c) 2018 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#include <numlib/math.h>
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

    SECTION("sign")
    {
        int x = 1;
        int y = -1;
        CHECK(sign(x, y) == -x);
    }

    SECTION("sqr")
    {
        int x = 2;
        CHECK(sqr(x) == 4);
    }
}
