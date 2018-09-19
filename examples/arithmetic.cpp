#include <iostream>
#include <numlib/matrix.h>

int main()
{
    using namespace Numlib;

    Mat<int> a = {{1, 2}, {3, 4}};
    Mat<int> b = {{10, 20}, {30, 40}};

    auto c = a + b;
    auto d = a * b;

    std::cout << "a = " << a << '\n'
              << "b = " << b << '\n'
              << "a + b = " << c << '\n'
              << "a * b = " << d << '\n';
}
