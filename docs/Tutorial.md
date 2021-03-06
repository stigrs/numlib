# Tutorial

## Table of Contents

* [Basic Matrix Uses](#basic-matrix-uses)
  + [Matrix Properties](#matrix-properties)
  + [Subscripting and Slicing](#subscripting-and-slicing)
  + [Sub-Matrix Views](#sub-matrix-views)
  + [Arithmetic Operations](#arithmetic-operations)
* [Constants](#constants)
* [Calculus](#calculus)
  + [Numerical Derivation](#numerical-derivation)
  + [Numerical Integration](#numerical-integration)
  + [Solve Ordinary Differential Equations](#solve-ordinary-differential-equations)
* [Linear Algebra](#linear-algebra)
  + [Basic Linear Algebra](#basic-linear-algebra)
  + [Matrix Inversion](#matrix-inversion)
  + [Matrix Decompositions](#matrix-decompositions)
  + [Eigensolvers](#eigensolvers)
  + [Solve Linear System of Equations](#solve-linear-system-of-equations)
  + [Linear Least Squares](#linear-least-squares)
  + [Band Matrices](#band-matrices)
  + [Sparse Matrices](#sparse-matrices)
    - [Eigensolver for Sparse Matrices](#eigensolver-for-sparse-matrices)
    - [Sparse Linear System of Equations](#sparse-linear-system-of-equations)
* [Special Functions](#special-functions)
* [Convolution](#convolution)
* [Discrete Fourier Transform](#discrete-fourier-transform)
* [Statistics](#statistics)

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

## Calculus

### Numerical Derivation
[back to top](#table-of-contents)

Example program:

    #include <iostream>
    #include <numlib/math.h>

    double f(double x) { return x * x; }

    int main() { std::cout << Numlib::dfdx(f, 2.0) << '\n'; }

Generated output:

    4

### Numerical Integration
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

        std::cout << "qags: " 
                  << qags([](double x, void*) { return std::sin(x); }, a, b) 
                  << '\n';
    }

Generated output:

    trapezoidal: 5.22
    5-point gaussian quadrature: 2
    qags: 2

### Solve Ordinary Differential Equations
[back to top](#table-of-contents)

Example program for nonstiff ODEs:

    #include <numlib/matrix.h>
    #include <numlib/math.h>
    #include <iostream>

    void lorenz(const Numlib::Vec<double>& x,
                Numlib::Vec<double>& dxdt,
                const double /* t */)
    {
        const double sigma(10.0);
        const double R(28.0);
        const double b(8.0 / 3.0);

        dxdt(0) = sigma * (x(1) - x(0));
        dxdt(1) = R * x(0) - x(1) - x(0) * x(2);
        dxdt(2) = -b * x(2) + x(0) * x(1);
    }

    int main()
    {
        Numlib::Vec<double> x = {10.0, 1.0, 1.0};

        double t0 = 0.0;
        double t1 = 0.1;

        for (int i = 0; i < 5; ++i) {
            Numlib::solve_ivp(lorenz, x, t0, t1);
            t1 += 0.1;
            std::cout << "At t = " << t0 << ", y = " << x(0) << " " 
                      << x(1) << " " << x(2) << '\n';
        }
    }

Generated output:

    At t = 0.1, y = 12.4201 22.1327 11.9965
    At t = 0.2, y = 19.5001 16.2247 45.2586
    At t = 0.3, y = 6.6136 -7.93158 37.7357
    At t = 0.4, y = -2.96399 -8.25055 28.2875
    At t = 0.5, y = -6.21703 -8.27847 25.1686

Example program for stiff ODEs:

    #include <numlib/math.h>
    #include <iostream>
    #include <iomanip>

    void fsys(double /* t */, double* y, double* ydot, void* /* data */)
    {
        ydot[0] = -0.04 * y[0] + 1.0e4 * y[1] * y[2];
        ydot[2] = 3.0e7 * y[1] * y[1];
        ydot[1] = -ydot[0] - ydot[2];
    }

    int main()
    {
        Numlib::Vec<double> rtol = {1.0e-4, 1.0e-8, 1.0e-4};
        Numlib::Vec<double> atol = {1.0e-6, 1.0e-10, 1.0e-6};

        Numlib::Odeint ode(fsys, rtol, atol);

        Numlib::Vec<double> y = {1.0, 0.0, 0.0};

        double t0 = 0.0;
        double t1 = 0.4;

        for (int i = 0; i < 12; ++i) {
            ode.integrate(y, t0, t1);
            std::cout << std::scientific << std::setprecision(6) 
                      << "At t = " << t0 << ", y = " << y(0) << " " 
                      << y(1) << " " << y(2) << '\n';
            t1 *= 10.0;
        }
    }

Generated output:

    At t = 4.000000e-01, y = 9.851712e-01 3.386380e-05 1.479493e-02
    At t = 4.000000e+00, y = 9.055333e-01 2.240655e-05 9.444430e-02
    At t = 4.000000e+01, y = 7.158403e-01 9.186334e-06 2.841505e-01
    At t = 4.000000e+02, y = 4.505250e-01 3.222964e-06 5.494717e-01
    At t = 4.000000e+03, y = 1.831976e-01 8.941773e-07 8.168015e-01
    At t = 4.000000e+04, y = 3.898729e-02 1.621940e-07 9.610125e-01
    At t = 4.000000e+05, y = 4.936362e-03 1.984221e-08 9.950636e-01
    At t = 4.000000e+06, y = 5.161833e-04 2.065787e-09 9.994838e-01
    At t = 4.000000e+07, y = 5.179804e-05 2.072027e-10 9.999482e-01
    At t = 4.000000e+08, y = 5.283686e-06 2.113485e-11 9.999947e-01
    At t = 4.000000e+09, y = 4.658646e-07 1.863459e-12 9.999995e-01
    At t = 4.000000e+10, y = 1.431682e-08 5.726733e-14 1.000000e+00

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

### Solve Linear System of Equations
[back to top](#table-of-contents)

Example program:

    #include <iostream>
    #include <numlib/matrix.h>
    #include <numlib/math.h>

    int main()
    {
        using namespace Numlib;

        Mat<double> A = {{1.0, 2.0, 3.0}, {2.0, 3.0, 4.0}, {3.0, 4.0, 1.0}};
        Mat<double> B = {{14.0}, {20.0}, {14.0}};

        linsolve(A, B);
        std::cout << B << '\n';
    }

Generated output:

    3 x 1
    [        1
             2
             3 ]

### Linear Least Squares
[back to top](#table-of-contents)

Example program:

    #include <iostream>
    #include <numlib/matrix.h>
    #include <numlib/math.h>

    int main()
    {
        using namespace Numlib;

        Mat<double> a = {{0.12, -8.19, 7.69, -2.26, -4.71},
                         {-6.91, 2.22, -5.12, -9.08, 9.96},
                         {-3.33, -8.94, -6.72, -4.40, -9.98},
                         {3.97, 3.33, -2.74, -7.92, -3.20}};

        Mat<double> b = {{7.30, 0.47, -6.28},
                         {1.33, 6.58, -3.42},
                         {2.68, -1.71, 3.46},
                         {-9.62, -0.79, 0.41},
                         {0.00, 0.00, 0.00}};

        lstsq(a, b);
        std::cout << b << '\n';
    }

Generated output:

    5 x 3
    [-0.685765 -0.243126 0.0642331
     -0.795005 -0.0836313  0.211787
      0.376681  0.120785 -0.654056
      0.288496 -0.241508  0.417589
      0.291609  0.352529 -0.301543 ]

### Band Matrices
[back to top](#table-of-contents)

Example program:

    #include <iostream>
    #include <numlib/matrix.h>
    #include <numlib/math.h>

    int main()
    {
        using namespace Numlib;

        Mat<double> a = {{1.0, 5.0, 2.0, 0.0},
                         {5.0, 2.0, 5.0, 2.0},
                         {2.0, 5.0, 3.0, 5.0},
                         {0.0, 2.0, 5.0, 4.0}};

        Band_matrix<double> ab(/* kl = */ 3, /* ku = */ 3, a);
        Mat<double> evec;
        Vec<double> eval;

        eigs(ab, evec, eval);
        std::cout << eval << '\n';
    }

Generated output:

    4
    [  -4.26201  -2.22987   3.95222   12.5397 ]

### Packed Matrices
[back to top](#table-of-contents)

Example program:

    #include <iostream>
    #include <numlib/matrix.h>
    #include <numlib/math.h>

    int main()
    {
    #ifdef USE_MKL
        using namespace Numlib;

        auto a = hilbert<>(5);
        Symm_mat<double, lo> ap(a);

        Mat<double> evec;
        Vec<double> eval(5); // correct size must be allocated

        eigs(ap, evec, eval);
        std::cout << eval << '\n';
    #else
        std::cout << "Intel MKL is required\n";
    #endif
    }

Generated output:

    5
    [ 3.28793e-06 0.000305898 0.0114075  0.208534   1.56705 ]

### Sparse Matrices

#### Eigensolver for Sparse Matrices
[back to top](#table-of-contents)

Example program:

    #include <iostream>
    #include <numlib/matrix.h>
    #include <numlib/math.h>

    int main()
    {
    #ifdef USE_MKL
        using namespace Numlib;

        // clang-format off
        BLAS_INT rows[12] = {0, 4, 9, 15, 22, 29, 36, 43, 50, 56, 61, 65};
        BLAS_INT cols[65] = {0,   1,   2,   3,
                             0,   1,   2,   3,   4,
                             0,   1,   2,   3,   4,   5,
                             0,   1,   2,   3,   4,   5,   6,
                                  1,   2,   3,   4,   5,   6,   7,
                                       2,   3,   4,   5,   6,   7,   8,
                                            3,   4,   5,   6,   7,   8,  9,
                                                 4,   5,   6,   7,   8,  9,  10,
                                                      5,   6,   7,   8,  9,  10,
                                                           6,   7,   8,  9,  10,
                                                                7,   8,  9,  10
        };
        double val[65] = {5.0, 2.0, 1.0, 1.0,
                          2.0, 6.0, 3.0, 1.0, 1.0,
                          1.0, 3.0, 6.0, 3.0, 1.0, 1.0,
                          1.0, 1.0, 3.0, 6.0, 3.0, 1.0, 1.0,
                               1.0, 1.0, 3.0, 6.0, 3.0, 1.0, 1.0,
                                    1.0, 1.0, 3.0, 6.0, 3.0, 1.0, 1.0,
                                         1.0, 1.0, 3.0, 6.0, 3.0, 1.0, 1.0,
                                              1.0, 1.0, 3.0, 6.0, 3.0, 1.0, 1.0,
                                                   1.0, 1.0, 3.0, 6.0, 3.0, 1.0,
                                                        1.0, 1.0, 3.0, 6.0, 2.0,
                                                             1.0, 1.0, 2.0, 5.0
        };
        // clang-format on

        Sp_mat<double> a(11, 11, val, cols, rows);
        Mat<double> evec(a.rows(), a.cols());
        Vec<double> eval(a.cols());

        double emin = 3.0;
        double emax = 7.0;

        eig(emin, emax, a, evec, eval);

        std::cout << eval << '\n';
    #else
        std::cout << "Intel MKL is required\n";
    #endif
    }

Generated output (in debug mode):

    Intel MKL Extended Eigensolvers: double precision driver
    Intel MKL Extended Eigensolvers: List of input parameters fpm(1:64)-- if different from default
    Intel MKL Extended Eigensolvers: fpm(1)=1
    Search interval [3.000000000000000e+00;7.000000000000000e+00]
    Intel MKL Extended Eigensolvers: Size subspace 11
    #Loop | #Eig  |    Trace     | Error-Trace |  Max-Residual
    Intel MKL Extended Eigensolvers: Resize subspace 9
    0,6,2.570747126011656e+01,1.000000000000000e+00,2.405853064980635e-07
    1,6,2.570747126011605e+01,7.308439567818174e-14,3.217497813885594e-14
    Intel MKL Extended Eigensolvers have successfully converged (to desired tolerance).
    Intel MKL Extended Eigensolvers have successfully converged (to desired tolerance).
    6
    [   3.17157         4         4   4.12925   4.40665         6 ]

#### Sparse Linear System of Equations
[back to top](#table-of-contents)

Example program:

    #include <iostream>
    #include <numlib/matrix.h>
    #include <numlib/math.h>

    int main()
    {
    #ifdef USE_MKL
        using namespace Numlib;

        Mat<double> A = {{0.0, 2.0, 0.0, 1.0, 0.0},
                         {4.0, -1.0, -1.0, 0.0, 0.0},
                         {0.0, 0.0, 0.0, 3.0, -6.0},
                         {-2.0, 0.0, 0.0, 0.0, 2.0},
                         {0.0, 0.0, 4.0, 2.0, 0.0}};

        Sp_mat<double> SA = gather(A);

        Mat<double> B = {{8.0}, {-1.0}, {-18.0}, {8.0}, {20.0}};
        Mat<double> x;

        linsolve(SA, B, x);
        std::cout << x << '\n';
    #else
        std::cout << "Intel MKL is required\n";
    #endif
    }

Generated output:

    5 x 1
    [        1
             2
             3
             4
             5 ]

## Special Functions
[back to top](#table-of-contents)

Example program:

    #include <iostream>
    #include <numlib/math.h>
    #include <complex>

    int main()
    {
        using namespace Faddeeva;

        std::complex<double> z(-1.0, 1.0);
        std::cout << "Faddeeva function: " << w(z) << '\n';

        z = {1.0, 2.0};
        std::cout << "Error function: " << erf(z) << '\n';

    #ifndef _MSC_VER
        z = {2.0, 1.0};
        std::cout << "Dawson function: " << Dawson(z) << '\n';
    #endif
    }

Generated output:

    Faddeeva function: (0.304744,-0.208219)
    Error function: (-0.536644,-5.04914)
    Dawson function: (0.163539,-0.153125)

## Convolution
[back to top](#table-of-contents)

Example program:

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

Generated output:

    6
    [         4        13        28        34        32        21 ]

## Discrete Fourier Transform
[back to top](#table-of-contents)

Example program:

    #include <iostream>
    #include <numlib/constants.h>
    #include <numlib/matrix.h>
    #include <numlib/math.h>

    double moda(int k, int l, int m)
    {
        return narrow_cast<double>(narrow_cast<long long>(k * l) % m);
    }

    int main()
    {
    #ifdef USE_MKL
        using namespace Numlib;

        const double twopi = 2.0 * Constants::pi;
        const int n = 6;
        int h = -1;

        Vec<double> x(n);

        double factor = (2 * (n - h) % n == 0) ? 1.0 : 2.0;
        for (int i = 0; i < n; ++i) {
            double phase = moda(i, h, n) / n;
            x(i) = factor * std::cos(twopi * phase) / n;
        }
        auto y = fft(x);
        std::cout << "x = " << x << "\n\n"
                  << "y = fft(x):\n"
                  << y << "\n\n"
                  << "x = ifft(y):\n"
                  << ifft(y) << '\n';

    #else
        std::cout << "Intel MKL is required\n";
    #endif
    }

Generated output:

    x = 6
    [  0.333333  0.166667 -0.166667 -0.333333 -0.166667  0.166667 ]
    
    y = fft(x):
    6
    [     (0,0) (1,-1.92296e-16) (0,1.92296e-16) (-2.22045e-16,0) (0,-1.92296e-16) (1,1.92296e-16) ]

    x = ifft(y):
    6
    [  0.333333  0.166667 -0.166667 -0.333333 -0.166667  0.166667 ]

## Statistics
[back to top](#table-of-contents)

Example program:

    #include <iostream>
    #include <numlib/matrix.h>
    #include <numlib/math.h>

    int main()
    {
        using namespace Numlib;

        Vec<double> a = {3.0,  13.0, 7.0,  5.0,  21.0, 23.0, 39.0, 23.0,
                         40.0, 23.0, 14.0, 12.0, 56.0, 23.0, 29.0};

        std::cout << "mean(a) =   " << mean(a) << '\n'
                  << "median(a) = " << median(a) << '\n'
                  << "stddev(a) = " << stddev(a) << '\n'
                  << "rms(a) =    " << rms(a) << "\n\n";

        Vec<double> b = {3.0,  13.0, 7.0,  5.0,  21.0, 23.0, 23.0,
                         40.0, 23.0, 14.0, 12.0, 56.0, 23.0, 29.0};
        Vec<double> c = {3.0,  13.0, 7.0,  5.0,  21.0, 23.0, 39.0,
                         23.0, 40.0, 23.0, 14.0, 12.0, 56.0, 23.0};

        std::cout << "cov(b, c) = " << cov(b, c) << "\n\n";

        Mat<double> p = {{0.8147, 0.9058, 0.1270}, {0.9134, 0.6324, 0.0975},
                         {0.2785, 0.5469, 0.9575}, {0.9649, 0.1576, 0.9706},
                         {0.9572, 0.4854, 0.8003}, {0.1419, 0.4218, 0.9157}};

        Mat<double> q = {{0.7922, 0.9595, 0.6557}, {0.0357, 0.8491, 0.9340},
                         {0.6787, 0.7577, 0.7431}, {0.3922, 0.6555, 0.1712},
                         {0.7060, 0.0318, 0.2769}, {0.0462, 0.0971, 0.8235}};

        std::cout << "rmsd(p, q) = " << kabsch_rmsd(p, q) << '\n';
    }

Generated output:

    mean(a) =   22.0667
    median(a) = 23
    stddev(a) = 14.4986
    rms(a) =    26.1368

    cov(b, c) = 59.7802

    rmsd(p, q) = 0.476097

[back to top](#table-of-contents)
