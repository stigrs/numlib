#include <iostream>
#include <numlib/matrix.h>
#include <numlib/math.h>

int main()
{
    using namespace Numlib;

    Vec<int> a = {1, 2, 3};
    Vec<int> b = {4, 5, 6, 7};

    std::cout << conv(a, b) << '\n';
    ;
}
