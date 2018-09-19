// Copyright (c) 2018 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#include <numlib/math.h>
#include <numlib/matrix.h>
#include <catch2/catch.hpp>

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
}
