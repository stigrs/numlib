# Tutorial

## Basic Matrix Uses

### Matrix Properties

Example program:

    #include <iostream>
    #include <numlib/matrix.h>

    int main()
    {
        using namespace Numlib;

        Mat<int> m1(3, 3);
        std::cout << "rank =      " << m1.rank() << '\n'
                  << "size =      " << m1.size() << '\n'
                  << "rows =      " << m1.rows() << '\n'
                  << "cols =      " << m1.cols() << '\n'
                  << "extent(0) = " << m1.extent(0) << '\n'
                  << "extent(1) = " << m1.extent(1) << "\n\n";

        m1.resize(5, 5);
        std::cout << "rank =      " << m1.rank() << '\n'
                  << "size =      " << m1.size() << '\n'
                  << "rows =      " << m1.rows() << '\n'
                  << "cols =      " << m1.cols() << '\n'
                  << "extent(0) = " << m1.extent(0) << '\n'
                  << "extent(1) = " << m1.extent(1) << "\n\n";

        Mat<int> m2 = {{1, 2, 3}, {4, 5, 6}};
        std::cout << m2 << '\n';
    }

Generated output:

    rank =      2
    size =      9
    rows =      3
    cols =      3
    extent(0) = 3
    extent(1) = 3

    rank =      2
    size =      25
    rows =      5
    cols =      5
    extent(0) = 5
    extent(1) = 5

    2 x 3
    [        1         2         3
             4         5         6 ]

### Subscripting and Slicing

Example program:

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

Generated output:

    index subscripting:
    v(0) = 1
    v(1) = 2
    v(2) = 3
    v(3) = 4
    v(4) = 5

    v(slice{0, 3}) (head):
    3
    [         1         2         3 ]

    v(slice{3}) (tail):
    2
    [         4         5 ]

    v(slice{0, 3, 2}) (slice with stride 2):
    3
    [         1         3         5 ]
