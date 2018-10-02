#include <iostream>
#include <numlib/matrix.h>
#include <numlib/math.h>

int main()
{
    using namespace Numlib;

    Mat<double> a = {{1.0, 5.0, 4.0, 2.0},
                     {-2.0, 3.0, 6.0, 4.0},
                     {5.0, 1.0, 0.0, -1.0},
                     {2.0, 3.0, -4.0, 0.0}};
    inv(a);
    std::cout << a << '\n';
}
