# Tutorial

## Table of Contents

* [Basic Matrix Uses](#basic-matrix-uses)
  + [Matrix Properties](#matrix-properties)
  + [Subscripting and Slicing](#subscripting-and-slicing)
  + [Sub-Matrix Views](#sub-matrix-views)
  + [Arithmetic Operations](#arithmetic-operations)
* [Constants](#constants)
* [Numerical Derivation](#numerical-derivation)
* [Numerical Integration](#numerical-integration)
* [Linear Algebra](#linear-algebra)
  + [Basic Linear Algebra](#basic-linear-algebra)
  + [Matrix Inversion](#matrix-inversion)
  + [Matrix Decompositions](#matrix-decompositions)
  + [Eigensolvers](#eigensolvers)

## Basic Matrix Uses

### Matrix Properties
[back to top](#table-of-contents)

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
        std::cout << "m2 = \n" << m2 << '\n';

        std::cout << "transpose(m2) =\n" << transpose(m2) << '\n';
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

    m2 =
    2 x 3
    [        1         2         3
             4         5         6 ]

    transpose(m2) =
    3 x 2
    [        1         4
             2         5
             3         6 ]

### Subscripting and Slicing
[back to top](#table-of-contents)

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
[back to top](#table-of-contents)

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
[back to top](#table-of-contents)

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
[back to top](#table-of-contents)

Example program:

    #include <iostream>
    #include <numlib/constants.h>

    int main() { std::cout << Numlib::Constants::pi << '\n'; }

Generated output:

    3.14159

## Numerical Derivation
[back to top](#table-of-contents)

Example program:

    #include <iostream>
    #include <numlib/math.h>

    double f(double x) { return x * x; }

    int main() { std::cout << Numlib::dfdx(f, 2.0) << '\n'; }

Generated output:

    4

## Numerical Integration
[back to top](#table-of-contents)

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

## Linear Algebra
 
### Basic Linear Algebra
[back to top](#table-of-contents)

Example program:

    #include <iostream>
    #include <numlib/matrix.h>
    #include <numlib/math.h>

    int main()
    {
        using namespace Numlib;

        Vec<double> x = linspace(1.0, 10.0, 10);
        std::cout << "x = \n" << x << "\n\n";

        Vec<double> y = 2.0 * ones<Vec<double>>(10);

        std::cout << "min(x) =    " << min(x) << '\n'
                  << "max(x) =    " << max(x) << '\n'
                  << "sum(x) =    " << sum(x) << '\n'
                  << "prod(x) =   " << prod(x) << '\n'
                  << "norm(x) =   " << norm(x) << '\n'
                  << "dot(x, y) = " << dot(x, y) << "\n\n";

        Vec<double> a = {3.0, -3.0, 1.0};
        Vec<double> b = {4.0, 9.0, 2.0};
        std::cout << "cross(a, b) = \n" << cross(a, b) << '\n';
    }

Generated output:

    x =
    10
    [         1         2         3         4         5         6         7
              8         9        10 ]

    min(x) =    1
    max(x) =    10
    sum(x) =    55
    prod(x) =   3.6288e+06
    norm(x) =   19.6214
    dot(x, y) = 110

    cross(a, b) =
    3
    [       -15        -2        39 ]

### Matrix Inversion
[back to top](#table-of-contents)

Example program:

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

Generated output:

    4 x 4
    [-0.190083  0.165289  0.280992 0.0578512
      0.347107 -0.214876 -0.165289 0.0247934
      0.165289 -0.0785124 0.0165289 -0.202479
     -0.603306   0.61157  0.239669   0.31405 ]

### Matrix Decompositions
[back to top](#table-of-contents)

Example program:

    #include <iostream>
    #include <numlib/matrix.h>
    #include <numlib/math.h>

    int main()
    {
        using namespace Numlib;

        Mat<double> a = {
            {12.0, -51.0, 4.0}, {6.0, 167.0, -68.0}, {-4.0, 24.0, -41.0}};

        Mat<double> q;
        Mat<double> r;

        qr(a, q, r);

        Mat<double> qr = q * r;
        std::cout << "QR decomposition (A = QR):\n" << qr << '\n';

        Mat<double> m = {
            {8.79, 9.93, 9.83, 5.45, 3.16},   
            {6.11, 6.91, 5.04, -0.27, 7.98},
            {-9.15, -7.93, 4.86, 4.85, 3.01}, 
            {9.57, 1.64, 8.83, 0.74, 5.8},
            {-3.49, 4.02, 9.80, 10.00, 4.27}, 
            {9.84, 0.15, -8.99, -6.02, -5.31}};

        Vec<double> s;
        Mat<double> u;
        Mat<double> vt;

        svd(m, s, u, vt);

        Mat<double> sigma = zeros<Mat<double>>(m.rows(), m.cols());
        sigma.diag() = s;

        std::cout << "SVD decomposition (M = USVT):\n" 
                  << u * sigma * vt 
                  << '\n';
    }

 Generated output:

    QR decomposition (A = QR):
    3 x 3
    [       12       -51         4
             6       167       -68
            -4        24       -41 ]

    SVD decomposition (M = USVT):
    6 x 5
    [     8.79      9.93      9.83      5.45      3.16
          6.11      6.91      5.04     -0.27      7.98
         -9.15     -7.93      4.86      4.85      3.01
          9.57      1.64      8.83      0.74       5.8
         -3.49      4.02       9.8        10      4.27
          9.84      0.15     -8.99     -6.02     -5.31 ]

### Eigensolvers
[back to top](#table-of-contents)

Example program:

    #include <iostream>
    #include <complex>
    #include <numlib/matrix.h>
    #include <numlib/math.h>

    int main()
    {
        using namespace Numlib;

        Mat<double> sa = hilbert(5);
        Vec<double> eval;

        eigs(sa, eval);
        std::cout << "Eigenvalues for symmetric matrix:\n"
                  << eval << "\n\n"
                  << "Eigenvectors for symmetric matrix:\n"
                  << sa << '\n';

        Mat<double> ga = {{1.0, 5.0, 4.0, 2.0},
                         {-2.0, 3.0, 6.0, 4.0},
                         {5.0, 1.0, 0.0, -1.0},
                         {2.0, 3.0, -4.0, 0.0}};

        Vec<std::complex<double>> geval;
        Mat<std::complex<double>> gevec;

        eig(ga, gevec, geval);
        std::cout << "Eigenvalues for general matrix:\n"
                  << geval << "\n\n"
                  << "Right eigenvectors for general matrix:\n"
                  << gevec << '\n';
    }

Generated output:

    Eigenvalues for symmetric matrix:
    5
    [ 3.28793e-06 0.000305898 0.0114075  0.208534   1.56705 ]

    Eigenvectors for symmetric matrix:
    5 x 5
    [-0.00617386 0.0471618  0.214214 -0.601871  0.767855
      0.116693 -0.432667 -0.724102  0.275913  0.445791
     -0.506164   0.66735 -0.120453  0.424877  0.321578
      0.767191  0.233025  0.309574  0.443903  0.253439
     -0.376246   -0.5576  0.565193  0.429013  0.209823 ]

    Eigenvalues for general matrix:
    4
    [ (-3.1736,1.12844) (-3.1736,-1.12844) (2.8422,0) (7.50501,0) ]

    Right eigenvectors for general matrix:
    4 x 4
    [(-0.168896,-0.112295) (-0.168896,0.112295) (-0.195144,0) (0.70846,0)
     (0.61502,-0.0394273) (0.61502,0.0394273) (0.0860169,0) (0.465904,0)
     (-0.19838,0.118805) (-0.19838,-0.118805) (-0.587648,0) (0.521106,0)
     (-0.724976,0) (-0.724976,-0) (0.780506,0) (0.0972959,0) ]
