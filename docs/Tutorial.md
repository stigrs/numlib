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

### Sub-Matrix Views

Example program:

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

Generated output:

    matrix:
    3 x 4
    [        1         2         3         4
             5         6         7         8
             9        10        11        12 ]

    row(1):
    4
    [         5         6         7         8 ]

    column(2):
    3
    [         3         7        11 ]

    m(0:2, 0:2):
    2 x 2
    [        1         2
             5         6 ]


    m(0:2:2, 1:2:2):
    2 x 2
    [        2         4
            10        12 ]


    matrix:
    3 x 4
    [        1         0         3         0
             5         6         7         8
             9         0        11         0 ]

### Arithmetic Operations

Example program:

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

Generated output:

    a = 2 x 2
    [        1         2
             3         4 ]

    b = 2 x 2
    [       10        20
            30        40 ]

    a + b = 2 x 2
    [       11        22
            33        44 ]

    a * b = 2 x 2
    [       70       100
           150       220 ]

## Constants

Example program:

    #include <iostream>
    #include <numlib/constants.h>

    int main() { std::cout << Numlib::Constants::pi << '\n'; }

Generated output:

    3.14159

## Numerical Derivation

Example program:

    #include <iostream>
    #include <numlib/math.h>

    double f(double x) { return x * x; }

    int main() { std::cout << Numlib::dfdx(f, 2.0) << '\n'; }

Generated output:

    4

## Numerical Integration

Example program:

    #include <iostream>
    #include <cmath>
    #include <numlib/matrix.h>
    #include <numlib/constants.h>
    #include <numlib/math.h>

    int main()
    {
        using namespace Numlib;

        double xlo = 2.1;
        double xup = 3.6;

        Vec<double> y = {3.2, 2.7, 2.9, 3.5, 4.1, 5.2};

        std::cout << "trapezoidal: " << trapz(xlo, xup, y) << '\n';

        double a = 0.0;
        double b = Constants::pi;

        std::cout << "5-point gaussian quadrature: "
                  << quad<5>([](double x) { return std::sin(x); }, a, b) 
                  << '\n';
    }

 Generated output:

    trapezoidal: 5.22
    5-point gaussian quadrature: 2
