// Copyright (c) 2018 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#include <numlib/matrix.h>
#include <catch2/catch.hpp>

TEST_CASE("test_sparse_matrix")
{
    using namespace Numlib;

    Mat<int> mat = {{1, 2, 0, 4, 0},
                    {6, 7, 0, 0, 0},
                    {0, 0, 13, 14, 15},
                    {16, 0, 18, 19, 0},
                    {0, 22, 0, 0, 25}};

    SECTION("element_access")
    {
        Sp_mat<int> spmat = gather(mat);
        CHECK(spmat(0, 0) == 1);
        CHECK(spmat(0, 2) == 0);
    }

    SECTION("insert")
    {
        Sp_mat<int> spmat = gather(mat);
        spmat.insert(0, 2, 3);
        spmat.insert(0, 4, 5);
        spmat.insert(1, 2, 8);
        spmat.insert(1, 3, 9);
        spmat.insert(1, 4, 10);
        spmat.insert(2, 0, 11);
        spmat.insert(2, 1, 12);
        spmat.insert(3, 1, 17);
        spmat.insert(3, 4, 20);
        spmat.insert(4, 0, 21);
        spmat.insert(4, 2, 23);
        spmat.insert(4, 3, 24);

        int iter = 1;
        for (int i = 0; i < 5; ++i) {
            for (int j = 0; j < 5; ++j) {
                CHECK(spmat(i, j) == iter);
                ++iter;
            }
        }
    }

    SECTION("scatter")
    {
        Sp_mat<int> spmat = gather(mat);
        CHECK(scatter(spmat) == mat);
    }

    SECTION("mv_mul")
    {
        Sp_mat<int> spmat = gather(mat);

        Vec<int> x = {1, 2, 3, 4, 5};
        Vec<int> ans = {21, 20, 170, 146, 169};

        CHECK(ans == spmat * x);
    }
}
