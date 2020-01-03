#include <iostream>
#include <numlib/matrix.h>

void beast(const Numlib::Matrix<Index, 2>& a)
{
    using namespace Numlib;

    // auto b = a(slice{5}, slice{5}); // does not compile
    Matrix<int, 2> b = a(slice{5}, slice{5});
    b(1, 1) = 666;
}

int main()
{
    using namespace Numlib;

    Matrix<Index, 2> a(5, 5);

    for (Index i = 0; i < a.rows(); ++i) {
        for (Index j = 0; j < a.cols(); ++j) {
            a(i, j) = i + j;
        }
    }
    std::cout << "Before:\n" << a << '\n';
    beast(a);
    std::cout << "After:\n" << a << '\n';
}
