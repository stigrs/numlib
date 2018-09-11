// Copyright (c) 2018 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#include <numlib/matrix.h>
#include <catch2/catch.hpp>

TEST_CASE("test_sparse_vector")
{
    using namespace Numlib;

    SECTION("element_access")
    {
        Sparse_vector<int> spvec = {{1, 10}, {4, 20}, {9, 30}};

        CHECK(spvec.num_nonzero() == 3);
        CHECK(spvec(1) == 10);
        CHECK(spvec(2) == 0);
        CHECK(spvec(3) == 0);
        CHECK(spvec(4) == 20);
        CHECK(spvec(5) == 0);
        CHECK(spvec(6) == 0);
        CHECK(spvec(7) == 0);
        CHECK(spvec(8) == 0);
        CHECK(spvec(9) == 30);
        CHECK(spvec(10) == 0);
    }

    SECTION("set_element")
    {
        Sparse_vector<int> spvec = {{1, 10}, {4, 20}, {9, 30}};
        spvec.values()[1] = 2;

        CHECK(spvec.num_nonzero() == 3);
        CHECK(spvec(1) == 10);
        CHECK(spvec(2) == 0);
        CHECK(spvec(3) == 0);
        CHECK(spvec(4) == 2);
        CHECK(spvec(5) == 0);
        CHECK(spvec(6) == 0);
        CHECK(spvec(7) == 0);
        CHECK(spvec(8) == 0);
        CHECK(spvec(9) == 30);
        CHECK(spvec(10) == 0);
    }

    SECTION("insert")
    {
        Sparse_vector<int> spvec = {{1, 10}, {4, 20}, {9, 30}};

        spvec.insert(5, 50);

        CHECK(spvec.num_nonzero() == 4);
        CHECK(spvec(1) == 10);
        CHECK(spvec(2) == 0);
        CHECK(spvec(3) == 0);
        CHECK(spvec(4) == 20);
        CHECK(spvec(5) == 50);
        CHECK(spvec(6) == 0);
        CHECK(spvec(7) == 0);
        CHECK(spvec(8) == 0);
        CHECK(spvec(9) == 30);
        CHECK(spvec(10) == 0);
    }
}
