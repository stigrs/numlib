// Copyright (c) 2018 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#ifndef NUMLIB_MATH_ODEINT_H
#define NUMLIB_MATH_ODEINT_H

#include <numlib/matrix.h>
#include <numlib/math_impl/lsoda.h>
#include <vector>

namespace Numlib {

// Class providing Livermore Solver for Ordinary Differential Equations
// with automatic method switching for stiff and nonstiff problems (LSODA).
//
// LSODA solves the initial value problem for stiff or nonstiff
// systems of first order ODEs:
//
//     dy/dt = f(t,y) ,  or, in component form,
//     dy(i)/dt = f(i) = f(i,t,y(0),y(1),...,y(neq-1)) (i = 0,...,neq-1).
//
// The function to be integrated must be of the form:
//
//     void f(double t, double* y, double* ydot, void* data);
//
class Odeint {
public:
    Odeint() = delete;

    Odeint(_lsoda_f fsys,
           double rtol_ = 1.0e-6,
           double atol_ = 1.0e-6,
           void* dptr = nullptr)
        : fptr(fsys), data_ptr(dptr)
    {
        rtol.push_back(0.0);
        atol.push_back(0.0);
        rtol.push_back(rtol_);
        atol.push_back(atol_);
        itol = 1;
        set_defaults();
    }

    Odeint(_lsoda_f fsys,
           const Vec<double>& rtol_,
           const Vec<double>& atol_,
           void* dptr = nullptr)
        : fptr(fsys), data_ptr(dptr)
    {
        rtol.push_back(0.0);
        atol.push_back(0.0);
        for (auto ri : rtol_) {
            rtol.push_back(ri);
        }
        for (auto ai : atol_) {
            atol.push_back(ai);
        }
        itol = 2;
        set_defaults();
    }

    // Destructor.
    ~Odeint() { n_lsoda_terminate(); }

    // Integrator.
    void integrate(Numlib::Vec<double>& y, double& t0, double& t1)
    {
        if (istate == 1) {
            set_init_value(y);
        }
        lsoda(fptr, neq, y0.data(), &t0, t1, itol, rtol.data(), atol.data(),
              itask, &istate, iopt, jt, iwork1, iwork2, iwork5, iwork6, iwork7,
              iwork8, iwork9, rwork1, rwork5, rwork6, rwork7, data_ptr);
        for (std::size_t i = 1; i < y0.size(); ++i) { // lsoda start index at 1
            y(i - 1) = y0[i];
        }
    }

private:
    void set_defaults();
    void set_init_value(const Vec<double>& yinit);

    void (*fptr)(double, double*, double*, void*);
    void* data_ptr;

    std::vector<double> rtol;
    std::vector<double> atol;
    std::vector<double> y0;

    int itol;
    int itask;
    int istate;
    int iopt;
    int jt;
    int neq;

    int iwork1;
    int iwork2;
    int iwork5;
    int iwork6;
    int iwork7;
    int iwork8;
    int iwork9;

    double rwork1;
    double rwork5;
    double rwork6;
    double rwork7;
};

inline void Odeint::set_defaults()
{
    itask = 1;
    istate = 1;
    iopt = 0;
    jt = 2;

    iwork1 = 0;
    iwork2 = 0;
    iwork5 = 0;
    iwork6 = 0;
    iwork7 = 0;
    iwork8 = 0;
    iwork9 = 0;

    rwork1 = 0.0;
    rwork5 = 0.0;
    rwork6 = 0.0;
    rwork7 = 0.0;
}

inline void Odeint::set_init_value(const Vec<double>& yinit)
{
    neq = narrow_cast<int>(yinit.size());
    y0.push_back(0.0); // note: lsoda start array index at 1
    for (auto yi : yinit) {
        y0.push_back(yi);
    }
}

} // namespace Numlib

#endif // NUMLIB_MATH_ODEINT_H
