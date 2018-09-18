#include <iostream>
#include <numlib/matrix.h>

int main()
{
    using namespace Numlib;

    Mat<int> m(3, 3);
    std::cout << "rank =      " << m.rank() << '\n'
              << "size =      " << m.size() << '\n'
              << "rows =      " << m.rows() << '\n'
              << "cols =      " << m.cols() << '\n'
              << "extent(0) = " << m.extent(0) << '\n'
              << "extent(1) = " << m.extent(1) << '\n';
}
