#include <iostream>
#include <numlib/matrix.h>

int main()
{
    using namespace Numlib;

    Matrix<int, 2> m = {{1, 2, 3, 4}, {5, 6, 7, 8}, {9, 10, 11, 12}};
    std::cout << "matrix:\n" << m << '\n';

    auto r = m.row(1);
    std::cout << "row(1):\n" << r << "\n\n";

    auto c = m.column(2);
    std::cout << "column(2):\n" << c << "\n\n";

    auto s1 = m(slice{0, 2}, slice{0, 2});
    std::cout << "m(0:2, 0:2):\n" << s1 << "\n\n";

    auto s2 = m(slice{0, 2, 2}, slice{1, 2, 2});
    std::cout << "m(0:2:2, 1:2:2):\n" << s2 << "\n\n";

    s2 = 0;
    std::cout << "matrix:\n" << m << '\n';
}
