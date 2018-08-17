////////////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2018 Stig Rune Sellevag. All rights reserved.
//
// This code is licensed under the MIT License (MIT).
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.
//
////////////////////////////////////////////////////////////////////////////////

#include <matrix/matrix.h>
#include <catch/catch.hpp>
#include <iostream>


TEST_CASE("test_matrix2")
{
    using namespace numlib;

    Matrix<int, 2> m2(3, 4);

    SECTION("size") { CHECK(m2.size() == 12); }

    SECTION("extents")
    {
        CHECK(m2.extent(0) == 3);
        CHECK(m2.extent(1) == 4);
    }

    SECTION("subscripting")
    {
        m2(0, 0) = 1;
        m2(0, 1) = 2;
        m2(0, 2) = 3;
        m2(0, 3) = 4;
        m2(1, 0) = 5;
        m2(1, 1) = 6;
        m2(1, 2) = 7;
        m2(1, 3) = 8;
        m2(2, 0) = 9;
        m2(2, 1) = 10;
        m2(2, 2) = 11;
        m2(2, 3) = 12;

        int it = 1;
        for (std::size_t i = 0; i < m2.extent(0); ++i) {
            for (std::size_t j = 0; j < m2.extent(1); ++j) {
                CHECK(m2(i, j) == it);
                ++it;
            }
        }
    }
}
