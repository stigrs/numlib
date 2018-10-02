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
