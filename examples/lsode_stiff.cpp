#include <numlib/math.h>
#include <iostream>

void my_fsys(int* /* neq */, double* /* t */, double* y, double* ydot)
{
    ydot[0] = -0.04 * y[0] + 1.0e4 * y[1] * y[2];
    ydot[2] = 3.0e7 * y[1] * y[1];
    ydot[1] = -ydot[0] - ydot[2];
}

void my_jsys(int* /* neq */,
             double* /* t */,
             double* y,
             int* /* ml */,
             int* /* mu */,
             double* pd,
             int* /* nrowpd */)
{
    // Column-major storage:
    pd[0] = -0.04;
    pd[1] = 0.04;
    pd[2] = 0.0;
    pd[3] = 1.0e4 * y[2];
    pd[5] = 6.0e7 * y[1];
    pd[4] = -pd[3] - pd[5];
    pd[6] = 1.0e4 * y[1];
    pd[7] = -pd[6];
    pd[8] = 0.0;
}

int main()
{
#ifdef ENABLE_LSODE
    Numlib::Lsode ode(my_fsys, my_jsys);

    Numlib::Vec<double> y = {1.0, 0.0, 0.0};

    double t0 = 0.0;
    double t1 = 0.4;

    for (int i = 0; i < 12; ++i) {
        ode.integrate(t0, t1, y);
        std::cout << "At t = " << t0 << ", y = " << y(0) << " " << y(1) << " "
                  << y(2) << '\n';
        t1 *= 10.0;
    }
#else
    std::cout << "Fortran compiler is needed\n";
#endif
}
