// Copyright (c) 2018 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#include <numlib/matrix.h>
#include <catch/catch.hpp>
#include <iostream>


TEST_CASE("test_matrix2")
{
    using namespace Numlib;

    Matrix<int, 2> m2 = {{1, 2, 3, 4}, {5, 6, 7, 8}, {9, 10, 11, 12}};

    SECTION("size") { CHECK(m2.size() == 12); }

    SECTION("extents")
    {
        CHECK(m2.extent(0) == 3);
        CHECK(m2.extent(1) == 4);
    }

    SECTION("subscripting")
    {
        int it = 1;
        for (std::size_t i = 0; i < m2.extent(0); ++i) {
            for (std::size_t j = 0; j < m2.extent(1); ++j) {
                CHECK(m2(i, j) == it);
                ++it;
            }
        }
    }
}
