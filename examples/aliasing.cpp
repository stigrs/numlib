#include <numlib/matrix.h>
#include <iostream>

using namespace Numlib;

int main()
{
    Mat<int> mat = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};

    std::cout << "Before:\n" << mat << '\n';
    mat(slice{1, 2}, slice{1, 2}) = mat(slice{0, 2}, slice{0, 2});
    std::cout << "After:\n" << mat << '\n';

    mat = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
    std::cout << "Before:\n" << mat << '\n';
    Mat<int> bv = mat(slice{0, 2}, slice{0, 2});
    mat(slice{1, 2}, slice{1, 2}) = bv;
    std::cout << "After:\n" << mat << '\n';
}
