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


TEST_CASE("test_matrix3")
{
    using namespace numlib;

    Matrix<int, 3> m3(2, 3, 4);

    SECTION("size") { CHECK(m3.size() == 24); }

    SECTION("extents")
    {
        CHECK(m3.extent(0) == 2);
        CHECK(m3.extent(1) == 3);
        CHECK(m3.extent(2) == 4);
    }

    SECTION("subscripting")
    {
        m3(0, 0, 0) = 1;
        m3(0, 0, 1) = 2;
        m3(0, 0, 2) = 3;
        m3(0, 0, 3) = 4;
        m3(0, 1, 0) = 5;
        m3(0, 1, 1) = 6;
        m3(0, 1, 2) = 7;
        m3(0, 1, 3) = 8;
        m3(0, 2, 0) = 9;
        m3(0, 2, 1) = 10;
        m3(0, 2, 2) = 11;
        m3(0, 2, 3) = 12;
        m3(1, 0, 0) = 13;
        m3(1, 0, 1) = 14;
        m3(1, 0, 2) = 15;
        m3(1, 0, 3) = 16;
        m3(1, 1, 0) = 17;
        m3(1, 1, 1) = 18;
        m3(1, 1, 2) = 19;
        m3(1, 1, 3) = 20;
        m3(1, 2, 0) = 21;
        m3(1, 2, 1) = 22;
        m3(1, 2, 2) = 23;
        m3(1, 2, 3) = 24;

        int it = 1;
        for (std::size_t i = 0; i < m3.extent(0); ++i) {
            for (std::size_t j = 0; j < m3.extent(1); ++j) {
                for (std::size_t k = 0; k < m3.extent(2); ++k) {
                    CHECK(m3(i, j, k) == it);
                    ++it;
                }
            }
        }
    }
}
