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

#include <numlib/matrix.h>
#include <catch/catch.hpp>
#include <iostream>


TEST_CASE("test_matrix3")
{
    using namespace Numlib;

    Matrix<int, 3> m3
        = {{{1, 2}, {3, 4}}, {{5, 6}, {7, 8}}, {{9, 10}, {11, 12}}};

    SECTION("size") { CHECK(m3.size() == 12); }

    SECTION("extents")
    {
        CHECK(m3.extent(0) == 3);
        CHECK(m3.extent(1) == 2);
        CHECK(m3.extent(2) == 2);
    }

    SECTION("subscripting")
    {
        CHECK(m3(0, 0, 0) == 1);
        CHECK(m3(0, 0, 1) == 2);
        CHECK(m3(0, 1, 0) == 3);
        CHECK(m3(0, 1, 1) == 4);
        CHECK(m3(1, 0, 0) == 5);
        CHECK(m3(1, 0, 1) == 6);
        CHECK(m3(1, 1, 0) == 7);
        CHECK(m3(1, 1, 1) == 8);
    }
}
