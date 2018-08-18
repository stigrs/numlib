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


TEST_CASE("test_matrix4")
{
    using namespace Numlib;

    Matrix<int, 4> m4(2, 3, 4, 5);

    SECTION("size") { CHECK(m4.size() == 120); }

    SECTION("extents")
    {
        CHECK(m4.extent(0) == 2);
        CHECK(m4.extent(1) == 3);
        CHECK(m4.extent(2) == 4);
        CHECK(m4.extent(3) == 5);
    }

    SECTION("subscripting")
    {
        m4(0, 0, 0, 0) = 1;
        m4(0, 0, 0, 1) = 2;

        CHECK(m4(0, 0, 0, 0) == 1);
        CHECK(m4(0, 0, 0, 1) == 2);
    }
}