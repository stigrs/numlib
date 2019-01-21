/*
  This is a C version of the LSODA library. I acquired the original
  source code from this web page:

    http://www.ccl.net/cca/software/SOURCES/C/kinetics2/index.shtml

  I merged several C files into one and added a simpler interface. I
  also made the array start from zero in functions called by lsoda(),
  and fixed two minor bugs: a) small memory leak in freevectors(); and
  b) misuse of lsoda() in the example.

  The original source code came with no license or copyright
  information. I now release this file under the MIT/X11 license. All
  authors' notes are kept in this file.

  - Heng Li <lh3lh3@gmail.com>
 */

/* The MIT License

   Copyright (c) 2009 Genome Research Ltd (GRL).

   Permission is hereby granted, free of charge, to any person obtaining
   a copy of this software and associated documentation files (the
   "Software"), to deal in the Software without restriction, including
   without limitation the rights to use, copy, modify, merge, publish,
   distribute, sublicense, and/or sell copies of the Software, and to
   permit persons to whom the Software is furnished to do so, subject to
   the following conditions:

   The above copyright notice and this permission notice shall be
   included in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.
*/

/* Contact: Heng Li <lh3@sanger.ac.uk> */

typedef void (*_lsoda_f)(double, double*, double*, void*);

/*
c-----------------------------------------------------------------------
c this is the march 30, 1987 version of
c lsoda.. livermore solver for ordinary differential equations, with
c         automatic method switching for stiff and nonstiff problems.
c
c this version is in double precision.
c
c lsoda solves the initial value problem for stiff or nonstiff
c systems of first order ode-s,
c     dy/dt = f(t,y) ,  or, in component form,
c     dy(i)/dt = f(i) = f(i,t,y(1),y(2),...,y(neq)) (i = 1,...,neq).
c
c this a variant version of the lsode package.
c it switches automatically between stiff and nonstiff methods.
c this means that the user does not have to determine whether the
c problem is stiff or not, and the solver will automatically choose the
c appropriate method.  it always starts with the nonstiff method.
c
c authors..
c                linda r. petzold  and  alan c. hindmarsh,
c                computing and mathematics research division, l-316
c                lawrence livermore national laboratory
c                livermore, ca 94550.
c
c references..
c 1.  alan c. hindmarsh,  odepack, a systematized collection of ode
c     solvers, in scientific computing, r. s. stepleman et al. (eds.),
c     north-holland, amsterdam, 1983, pp. 55-64.
c 2.  linda r. petzold, automatic selection of methods for solving
c     stiff and nonstiff systems of ordinary differential equations,
c     siam j. sci. stat. comput. 4 (1983), pp. 136-148.
c-----------------------------------------------------------------------
c summary of usage.
c
c communication between the user and the lsoda package, for normal
c situations, is summarized here.  this summary describes only a subset
c of the full set of options available.  see the full description for
c details, including alternative treatment of the jacobian matrix,
c optional inputs and outputs, nonstandard options, and
c instructions for special situations.  see also the example
c problem (with program and output) following this summary.
c
c a. first provide a subroutine of the form..
c               subroutine f (neq, t, y, ydot)
c               dimension y(neq), ydot(neq)
c which supplies the vector function f by loading ydot(i) with f(i).
c
c b. write a main program which calls subroutine lsoda once for
c each point at which answers are desired.  this should also provide
c for possible use of logical unit 6 for output of error messages
c by lsoda.  on the first call to lsoda, supply arguments as follows..
c f      = name of subroutine for right-hand side vector f.
c          this name must be declared external in calling program.
c neq    = number of first order ode-s.
c y      = array of initial values, of length neq.
c t      = the initial value of the independent variable.
c tout   = first point where output is desired (.ne. t).
c itol   = 1 or 2 according as atol (below) is a scalar or array.
c rtol   = relative tolerance parameter (scalar).
c atol   = absolute tolerance parameter (scalar or array).
c          the estimated local error in y(i) will be controlled so as
c          to be less than
c             ewt(i) = rtol*abs(y(i)) + atol     if itol = 1, or
c             ewt(i) = rtol*abs(y(i)) + atol(i)  if itol = 2.
c          thus the local error test passes if, in each component,
c          either the absolute error is less than atol (or atol(i)),
c          or the relative error is less than rtol.
c          use rtol = 0.0 for pure absolute error control, and
c          use atol = 0.0 (or atol(i) = 0.0) for pure relative error
c          control.  caution.. actual (global) errors may exceed these
c          local tolerances, so choose them conservatively.
c itask  = 1 for normal computation of output values of y at t = tout.
c istate = integer flag (input and output).  set istate = 1.
c iopt   = 0 to indicate no optional inputs used.
c rwork  = real work array of length at least..
c             22 + neq * max(16, neq + 9).
c          see also paragraph e below.
c lrw    = declared length of rwork (in user-s dimension).
c iwork  = integer work array of length at least  20 + neq.
c liw    = declared length of iwork (in user-s dimension).
c jac    = name of subroutine for jacobian matrix.
c          use a dummy name.  see also paragraph e below.
c jt     = jacobian type indicator.  set jt = 2.
c          see also paragraph e below.
c note that the main program must declare arrays y, rwork, iwork,
c and possibly atol.
c
c c. the output from the first call (or any call) is..
c      y = array of computed values of y(t) vector.
c      t = corresponding value of independent variable (normally tout).
c istate = 2  if lsoda was successful, negative otherwise.
c          -1 means excess work done on this call (perhaps wrong jt).
c          -2 means excess accuracy requested (tolerances too small).
c          -3 means illegal input detected (see printed message).
c          -4 means repeated error test failures (check all inputs).
c          -5 means repeated convergence failures (perhaps bad jacobian
c             supplied or wrong choice of jt or tolerances).
c          -6 means error weight became zero during problem. (solution
c             component i vanished, and atol or atol(i) = 0.)
c          -7 means work space insufficient to finish (see messages).
c
c d. to continue the integration after a successful return, simply
c reset tout and call lsoda again.  no other parameters need be reset.
c
c e. note.. if and when lsoda regards the problem as stiff, and
c switches methods accordingly, it must make use of the neq by neq
c jacobian matrix, j = df/dy.  for the sake of simplicity, the
c inputs to lsoda recommended in paragraph b above cause lsoda to
c treat j as a full matrix, and to approximate it internally by
c difference quotients.  alternatively, j can be treated as a band
c matrix (with great potential reduction in the size of the rwork
c array).  also, in either the full or banded case, the user can supply
c j in closed form, with a routine whose name is passed as the jac
c argument.  these alternatives are described in the paragraphs on
c rwork, jac, and jt in the full description of the call sequence below.
c
c-----------------------------------------------------------------------
*/

#ifndef NUMLIB_MATH_LSODA_H
#define NUMLIB_MATH_LSODA_H

#ifdef __cplusplus
extern "C" {
#endif

void lsoda(_lsoda_f f,
           int neq,
           double* y,
           double* t,
           double tout,
           int itol,
           double* rtol,
           double* atol,
           int itask,
           int* istate,
           int iopt,
           int jt,
           int iwork1,
           int iwork2,
           int iwork5,
           int iwork6,
           int iwork7,
           int iwork8,
           int iwork9,
           double rwork1,
           double rwork5,
           double rwork6,
           double rwork7,
           void* _data);

void n_lsoda_terminate(void);

#ifdef __cplusplus
}
#endif

#endif /* NUMLIB_MATH_LSODA_H */
