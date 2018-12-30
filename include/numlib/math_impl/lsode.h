// Copyright (c) 2018 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#ifndef NUMLIB_MATH_LSODE_H
#define NUMLIB_MATH_LSODE_H

#include <numlib/math_impl/dlsode.h>

#ifdef ENABLE_ODESOLVERS

namespace Numlib {

// Enumeration of Lsode method flags.
enum Lsode_methods { nonstiff = 10, stiff_user_jac = 21, stiff_int_jac = 22 };

// Class providing Livermore Solver for Ordinary Differential Equations.
//
// Note:
// -----
// Currently,
class Lsode {
public:
    // Constructor.
    Lsode(void (*fsys)(int* neq, double* t, double* y, double* ydot),
          void (*jsys)(int* neq,
                       double* t,
                       double* y,
                       int* ml,
                       int* mu,
                       double* pd,
                       int* nrowpd),
          int n,
          Lsode_methods flag_ = nonstiff,
          double rtol_ = 1.0e-6,
          double atol_ = 1.0e-6)
        : fptr(fsys), jptr(jsys), neq{n}, flag{flag_}, rtol{rtol_}, atol{atol_}
    {
        set_defaults();
        allocate_memory();
    }

    // Destructor.
    ~Lsode()
    {
        if (lrw > 0) {
            delete[] rwork;
        }
        if (liw > 0) {
            delete[] iwork;
        }
    }

    // Integrator.
    void integrate(double& t0, double& t1, Numlib::Vec<double>& y)
    {
        yptr = y.data();

        dlsode_(fptr, &neq, yptr, &t0, &t1, &itol, &rtol, &atol, &itask,
                &istate, &iopt, rwork, &lrw, iwork, &liw, jptr, &mf);
    }

private:
    void set_defaults();
    void allocate_memory();

    void (*fptr)(int* neq, double* t, double* y, double* ydot);
    void (*jptr)(int* neq,
                 double* t,
                 double* y,
                 int* ml,
                 int* mu,
                 double* pd,
                 int* nrowpd);

    int neq;

    Lsode_methods flag;

    double rtol;
    double atol;

    int mf;
    int itask;
    int istate;
    int itol;
    int iopt;

    int lrw;
    int liw;

    double* yptr;
    double* rwork;
    int* iwork;
};

inline void Lsode::set_defaults()
{
    itask = 1;
    istate = 1;
    itol = 1;
    iopt = 0;
}

inline void Lsode::allocate_memory()
{
    switch (flag) {
    case stiff_user_jac:
        mf = 21;
        lrw = 22 + 9 * neq + neq * neq;
        liw = 20 + neq;
        break;
    case stiff_int_jac:
        mf = 22;
        lrw = 22 + 9 * neq + neq * neq;
        liw = 20 + neq;
        break;
    case nonstiff:
    default:
        mf = 10;
        lrw = 20 + 16 * neq;
        liw = 20;
    }
    rwork = new double[lrw];
    iwork = new int[liw];
}

} // namespace Numlib

#endif // ENABLE_ODESOLVERS

#endif // NUMLIB_MATH_LSODE_H
