// Copyright (c) 2018 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#include <numlib/matrix.h>
#include <catch2/catch.hpp>

TEST_CASE("test_packed_matrix")
{
    using namespace Numlib;

    SECTION("upper_triangangular")
    {
        // | 1 4 1 |
        // | 0 6 4 |
        // | 0 0 1 |

        int upper[6] = {1, 4, 1, 6, 4, 1};
        Upper_triang_mat<int> u(3, upper);

        CHECK(u.uplo_scheme() == 'U');

        CHECK(u(0, 0) == 1);
        CHECK(u(0, 1) == 4);
        CHECK(u(0, 2) == 1);
        CHECK(u(1, 1) == 6);
        CHECK(u(1, 2) == 4);
        CHECK(u(2, 2) == 1);
    }

    SECTION("lower_triangangular")
    {
        // | 1 0 0 |
        // | 2 8 0 |
        // | 4 9 7 |

        Mat<int> lower = {{1, 0, 0}, {2, 8, 0}, {4, 9, 7}};
        Lower_triang_mat<int> l(lower);

        CHECK(l.uplo_scheme() == 'L');

        CHECK(l(0, 0) == 1);
        CHECK(l(1, 0) == 2);
        CHECK(l(1, 1) == 8);
        CHECK(l(2, 0) == 4);
        CHECK(l(2, 1) == 9);
        CHECK(l(2, 2) == 7);
    }

    SECTION("swap")
    {
        // | 1 0 0 |
        // | 2 8 0 |
        // | 4 9 7 |

        int l1[6] = {1, 2, 8, 4, 9, 7};
        Packed_matrix<int> la(3, l1);

        // | 1 0 0 |
        // | 4 1 0 |
        // | 6 4 1 |

        int l2[6] = {1, 4, 1, 6, 4, 1};
        Packed_matrix<int> lb(3, l2);

        std::swap(la, lb);

        CHECK(la(0, 0) == 1);
        CHECK(la(1, 0) == 4);
        CHECK(la(1, 1) == 1);
        CHECK(la(2, 0) == 6);
        CHECK(la(2, 1) == 4);
        CHECK(la(2, 2) == 1);
    }

    SECTION("operator+=")
    {
        // | 1 4 1 |
        // | 0 6 4 |
        // | 0 0 1 |

        Mat<int> a = {{1, 4, 1}, {0, 6, 4}, {0, 0, 1}};

        Upper_triang_mat<int> u(a);

        u += 1;

        CHECK(u(0, 0) == 2);
        CHECK(u(0, 1) == 5);
        CHECK(u(0, 2) == 2);
        CHECK(u(1, 1) == 7);
        CHECK(u(1, 2) == 5);
        CHECK(u(2, 2) == 2);
    }
}
