#include <iostream>
#include <numlib/matrix.h>

int main()
{
    using namespace Numlib;

    Vec<int> v = {1, 2, 3, 4, 5};

    std::cout << "index subscripting:\n";
    for (Index i = 0; i < v.size(); ++i) {
        std::cout << "v(" << i << ") = " << v(i) << '\n';
    }
    std::cout << '\n';

    std::cout << "v(slice{0, 3}) (head):\n"
              << v(slice{0, 3}) << "\n\n"
              << "v(slice{3}) (tail):\n"
              << v(slice{3}) << "\n\n"
              << "v(slice{0, 3, 2}) (slice with stride 2):\n"
              << v(slice{0, 3, 2}) << '\n';
}
