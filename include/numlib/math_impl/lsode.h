// Copyright (c) 2018 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#ifndef NUMLIB_MATH_LSODE_H
#define NUMLIB_MATH_LSODE_H

#ifdef ENABLE_ODEPACK
#include <numlib/math_impl/odepack.h>

namespace Numlib {

// Enumeration of Lsode method flags.
enum Lsode_method {
    ode_nonstiff = 10,
    ode_stiff_user_jac = 21,
    ode_stiff_int_jac = 22
};

// Class providing Livermore Solver for Ordinary Differential Equations.
//
// Note:
// -----
// Currently, only default options are supported as well as dense
// Jacobians.
//
class Lsode {
public:
    Lsode() = delete;

    // Constructor for nonstiff ODEs or stiff ODEs with generated Jacobian.
    Lsode(void (*fsys)(int& neq, double& t, double* y, double* ydot),
          Lsode_method flag_ = ode_nonstiff,
          double rtol_ = 1.0e-6,
          double atol_ = 1.0e-6)
        : fptr(fsys), flag{flag_}, rtol{rtol_}, atol{atol_}
    {
        jptr = jdummy;
        flag = ode_stiff_user_jac;
        set_defaults();
    }

    // Constructor for stiff ODEs with user-supplied Jacobian.
    Lsode(void (*fsys)(int& neq, double& t, double* y, double* ydot),
          void (*jsys)(int& neq,
                       double& t,
                       double* y,
                       int& ml,
                       int& mu,
                       double* pd,
                       int& nrowpd),
          double rtol_ = 1.0e-6,
          double atol_ = 1.0e-6)
        : fptr(fsys), jptr(jsys), rtol{rtol_}, atol{atol_}
    {
        flag = ode_stiff_user_jac;
        set_defaults();
    }

    // Destructor.
    ~Lsode() = default;

    // Integrator.
    void integrate(Numlib::Vec<double>& y, double& t0, double& t1)
    {
        if (istate == 1) {
            allocate_memory(narrow_cast<int>(y.size()));
        }
        dlsode_(fptr, neq, y.data(), t0, t1, itol, &rtol, &atol, itask, istate,
                iopt, rwork.data(), lrw, iwork.data(), liw, jptr, mf);
    }

private:
    static void jdummy(int& /* neq */,
                       double& /* t */,
                       double* /* y */,
                       int& /* ml */,
                       int& /* mu */,
                       double* /* pd */,
                       int& /* nrowpd */)
    {
        // do nothing
    }

    void set_defaults();
    void allocate_memory(int neq_);

    void (*fptr)(int& neq, double& t, double* y, double* ydot);
    void (*jptr)(int& neq,
                 double& t,
                 double* y,
                 int& ml,
                 int& mu,
                 double* pd,
                 int& nrowpd);

    Lsode_method flag;

    double rtol;
    double atol;

    int mf;
    int itask;
    int istate;
    int itol;
    int iopt;

    int neq;
    int lrw;
    int liw;

    Vec<double> rwork;
    Vec<int> iwork;
};

inline void Lsode::set_defaults()
{
    itask = 1;
    istate = 1;
    itol = 1;
    iopt = 0;
}

inline void Lsode::allocate_memory(int neq_)
{
    neq = neq_;

    switch (flag) {
    case ode_stiff_user_jac:
        mf = 21;
        lrw = 22 + 9 * neq + neq * neq;
        liw = 20 + neq;
        break;
    case ode_stiff_int_jac:
        mf = 22;
        lrw = 22 + 9 * neq + neq * neq;
        liw = 20 + neq;
        break;
    case ode_nonstiff:
    default:
        mf = 10;
        lrw = 20 + 16 * neq;
        liw = 20;
    }
    rwork.resize(lrw);
    iwork.resize(liw);
}

} // namespace Numlib

#endif // ENABLE_ODEPACK

#endif // NUMLIB_MATH_LSODE_H
