      subroutine rkf45(f,neqn,y,t,tout,relerr,abserr,iflag,work,iwork)
c
c     fehlberg fourth-fifth order runge-kutta method
c
c     written by h.a.watts and l.f.shampine
c                   sandia laboratories
c                  albuquerque,new mexico
c
c    rkf45 is primarily designed to solve non-stiff and mildly stiff
c    differential equations when derivative evaluations are inexpensive.
c    rkf45 should generally not be used when the user is demanding
c    high accuracy.
c
c abstract
c
c    subroutine  rkf45  integrates a system of neqn first order
c    ordinary differential equations of the form
c             dy(i)/dt = f(t,y(1),y(2),...,y(neqn))
c              where the y(i) are given at t .
c    typically the subroutine is used to integrate from t to tout but it
c    can be used as a one-step integrator to advance the solution a
c    single step in the direction of tout.  on return the parameters in
c    the call list are set for continuing the integration. the user has
c    only to call rkf45 again (and perhaps define a new value for tout).
c    actually, rkf45 is an interfacing routine which calls subroutine
c    rkfs for the solution.  rkfs in turn calls subroutine  fehl which
c    computes an approximate solution over one step.
c
c    rkf45  uses the runge-kutta-fehlberg (4,5)  method described
c    in the reference
c    e.fehlberg , low-order classical runge-kutta formulas with stepsize
c                 control , nasa tr r-315
c
c    the performance of rkf45 is illustrated in the reference
c    l.f.shampine,h.a.watts,s.davenport, solving non-stiff ordinary
c                 differential equations-the state of the art ,
c                 sandia laboratories report sand75-0182 ,
c                 to appear in siam review.
c
c
c    the parameters represent-
c      f -- subroutine f(t,y,yp) to evaluate derivatives yp(i)=dy(i)/dt
c      neqn -- number of equations to be integrated
c      y(*) -- solution vector at t
c      t -- independent variable
c      tout -- output point at which solution is desired
c      relerr,abserr -- relative and absolute error tolerances for local
c            error test. at each step the code requires that
c                 abs(local error) .le. relerr*abs(y) + abserr
c            for each component of the local error and solution vectors
c      iflag -- indicator for status of integration
c      work(*) -- array to hold information internal to rkf45 which is
c            necessary for subsequent calls. must be dimensioned
c            at least  3+6*neqn
c      iwork(*) -- integer array used to hold information internal to
c            rkf45 which is necessary for subsequent calls. must be
c            dimensioned at least  5
c
c
c  first call to rkf45
c
c    the user must provide storage in his calling program for the arrays
c    in the call list  -      y(neqn) , work(3+6*neqn) , iwork(5)  ,
c    declare f in an external statement, supply subroutine f(t,y,yp) and
c    initialize the following parameters-
c
c      neqn -- number of equations to be integrated.  (neqn .ge. 1)
c      y(*) -- vector of initial conditions
c      t -- starting point of integration , must be a variable
c      tout -- output point at which solution is desired.
c            t=tout is allowed on the first call only, in which case
c            rkf45 returns with iflag=2 if continuation is possible.
c      relerr,abserr -- relative and absolute local error tolerances
c            which must be non-negative. relerr must be a variable while
c            abserr may be a constant. the code should normally not be
c            used with relative error control smaller than about 1.e-8 .
c            to avoid limiting precision difficulties the code requires
c            relerr to be larger than an internally computed relative
c            error parameter which is machine dependent. in particular,
c            pure absolute error is not permitted. if a smaller than
c            allowable value of relerr is attempted, rkf45 increases
c            relerr appropriately and returns control to the user before
c            continuing the integration.
c      iflag -- +1,-1  indicator to initialize the code for each new
c            problem. normal input is +1. the user should set iflag=-1
c            only when one-step integrator control is essential. in this
c            case, rkf45 attempts to advance the solution a single step
c            in the direction of tout each time it is called. since this
c            mode of operation results in extra computing overhead, it
c            should be avoided unless needed.
c
c
c  output from rkf45
c
c      y(*) -- solution at t
c      t -- last point reached in integration.
c      iflag = 2 -- integration reached tout. indicates successful retur
c                   and is the normal mode for continuing integration.
c            =-2 -- a single successful step in the direction of tout
c                   has been taken. normal mode for continuing
c                   integration one step at a time.
c            = 3 -- integration was not completed because relative error
c                   tolerance was too small. relerr has been increased
c                   appropriately for continuing.
c            = 4 -- integration was not completed because more than
c                   3000 derivative evaluations were needed. this
c                   is approximately 500 steps.
c            = 5 -- integration was not completed because solution
c                   vanished making a pure relative error test
c                   impossible. must use non-zero abserr to continue.
c                   using the one-step integration mode for one step
c                   is a good way to proceed.
c            = 6 -- integration was not completed because requested
c                   accuracy could not be achieved using smallest
c                   allowable stepsize. user must increase the error
c                   tolerance before continued integration can be
c                   attempted.
c            = 7 -- it is likely that rkf45 is inefficient for solving
c                   this problem. too much output is restricting the
c                   natural stepsize choice. use the one-step integrator
c                   mode.
c            = 8 -- invalid input parameters
c                   this indicator occurs if any of the following is
c                   satisfied -   neqn .le. 0
c                                 t=tout  and  iflag .ne. +1 or -1
c                                 relerr or abserr .lt. 0.
c                                 iflag .eq. 0  or  .lt. -2  or  .gt. 8
c      work(*),iwork(*) -- information which is usually of no interest
c                   to the user but necessary for subsequent calls.
c                   work(1),...,work(neqn) contain the first derivatives
c                   of the solution vector y at t. work(neqn+1) contains
c                   the stepsize h to be attempted on the next step.
c                   iwork(1) contains the derivative evaluation counter.
c
c
c  subsequent calls to rkf45
c
c    subroutine rkf45 returns with all information needed to continue
c    the integration. if the integration reached tout, the user need onl
c    define a new tout and call rkf45 again. in the one-step integrator
c    mode (iflag=-2) the user must keep in mind that each step taken is
c    in the direction of the current tout. upon reaching tout (indicated
c    by changing iflag to 2),the user must then define a new tout and
c    reset iflag to -2 to continue in the one-step integrator mode.
c
c    if the integration was not completed but the user still wants to
c    continue (iflag=3,4 cases), he just calls rkf45 again. with iflag=3
c    the relerr parameter has been adjusted appropriately for continuing
c    the integration. in the case of iflag=4 the function counter will
c    be reset to 0 and another 3000 function evaluations are allowed.
c
c    however,in the case iflag=5, the user must first alter the error
c    criterion to use a positive value of abserr before integration can
c    proceed. if he does not,execution is terminated.
c
c    also,in the case iflag=6, it is necessary for the user to reset
c    iflag to 2 (or -2 when the one-step integration mode is being used)
c    as well as increasing either abserr,relerr or both before the
c    integration can be continued. if this is not done, execution will
c    be terminated. the occurrence of iflag=6 indicates a trouble spot
c    (solution is changing rapidly,singularity may be present) and it
c    often is inadvisable to continue.
c
c    if iflag=7 is encountered, the user should use the one-step
c    integration mode with the stepsize determined by the code or
c    consider switching to the adams codes de/step,intrp. if the user
c    insists upon continuing the integration with rkf45, he must reset
c    iflag to 2 before calling rkf45 again. otherwise,execution will be
c    terminated.
c
c    if iflag=8 is obtained, integration can not be continued unless
c    the invalid input parameters are corrected.
c
c    it should be noted that the arrays work,iwork contain information
c    required for subsequent integration. accordingly, work and iwork
c    should not be altered.
c
c
      integer neqn,iflag,iwork(5)
      double precision y(neqn),t,tout,relerr,abserr,work(1)
c
      external f
c
      integer k1,k2,k3,k4,k5,k6,k1m
c
c
c     compute indices for the splitting of the work array
c
      k1m=neqn+1
      k1=k1m+1
      k2=k1+neqn
      k3=k2+neqn
      k4=k3+neqn
      k5=k4+neqn
      k6=k5+neqn
c
c     this interfacing routine merely relieves the user of a long
c     calling list via the splitting apart of two working storage
c     arrays. if this is not compatible with the users compiler,
c     he must use rkfs directly.
c
      call rkfs(f,neqn,y,t,tout,relerr,abserr,iflag,work(1),work(k1m),
     1          work(k1),work(k2),work(k3),work(k4),work(k5),work(k6),
     2          work(k6+1),iwork(1),iwork(2),iwork(3),iwork(4),iwork(5))
c
      return
      end
      subroutine rkfs(f,neqn,y,t,tout,relerr,abserr,iflag,yp,h,f1,f2,f3,
     1                f4,f5,savre,savae,nfe,kop,init,jflag,kflag)
c
c     fehlberg fourth-fifth order runge-kutta method
c
c
c     rkfs integrates a system of first order ordinary differential
c     equations as described in the comments for rkf45 .
c     the arrays yp,f1,f2,f3,f4,and f5 (of dimension at least neqn) and
c     the variables h,savre,savae,nfe,kop,init,jflag,and kflag are used
c     internally by the code and appear in the call list to eliminate
c     local retention of variables between calls. accordingly, they
c     should not be altered. items of possible interest are
c         yp - derivative of solution vector at t
c         h  - an appropriate stepsize to be used for the next step
c         nfe- counter on the number of derivative function evaluations
c
c
      logical hfaild,output
c
      integer  neqn,iflag,nfe,kop,init,jflag,kflag
      double precision  y(neqn),t,tout,relerr,abserr,h,yp(neqn),
     1  f1(neqn),f2(neqn),f3(neqn),f4(neqn),f5(neqn),savre,
     2  savae
c
      external f
c
      double precision  a,ae,dt,ee,eeoet,esttol,et,hmin,remin,rer,s,
     1  scale,tol,toln,twoeps,u26,ypk
c
      integer  k,maxnfe,mflag
c
      double precision  dabs,dmax1,dmin1,dsign,d1mach
c
c  remin is the minimum acceptable value of relerr.  attempts
c  to obtain higher accuracy with this subroutine are usually
c  very expensive and often unsuccessful.
c
      data remin/1.d-12/
c
c
c     the expense is controlled by restricting the number
c     of function evaluations to be approximately maxnfe.
c     as set, this corresponds to about 500 steps.
c
      data maxnfe/3000/
c
c   here two constants emboding the machine epsilon is present
c   twoesp is set to twice the machine epsilon while u26 is set
c   to 26 times the machine epsilon
c
c     data twoeps, u26/4.4d-16, 5.72d-15/                               ***
      twoeps = 2.*d1mach(4)                                             ***
      u26 = 13.*twoeps                                                  ***
c
c
c     check input parameters
c
c
      if (neqn .lt. 1) go to 10
      if ((relerr .lt. 0.0d0)  .or.  (abserr .lt. 0.0d0)) go to 10
      mflag=iabs(iflag)
      if ((mflag .ge. 1)  .and.  (mflag .le. 8)) go to 20
c
c     invalid input
   10 iflag=8
      return
c
c     is this the first call
   20 if (mflag .eq. 1) go to 50
c
c     check continuation possibilities
c
      if ((t .eq. tout) .and. (kflag .ne. 3)) go to 10
      if (mflag .ne. 2) go to 25
c
c     iflag = +2 or -2
      if (kflag .eq. 3) go to 45
      if (init .eq. 0) go to 45
      if (kflag .eq. 4) go to 40
      if ((kflag .eq. 5)  .and.  (abserr .eq. 0.0d0)) go to 30
      if ((kflag .eq. 6)  .and.  (relerr .le. savre)  .and.
     1    (abserr .le. savae)) go to 30
      go to 50
c
c     iflag = 3,4,5,6,7 or 8
   25 if (iflag .eq. 3) go to 45
      if (iflag .eq. 4) go to 40
      if ((iflag .eq. 5) .and. (abserr .gt. 0.0d0)) go to 45
c
c     integration cannot be continued since user did not respond to
c     the instructions pertaining to iflag=5,6,7 or 8
   30 stop
c
c     reset function evaluation counter
   40 nfe=0
      if (mflag .eq. 2) go to 50
c
c     reset flag value from previous call
   45 iflag=jflag
      if (kflag .eq. 3) mflag=iabs(iflag)
c
c     save input iflag and set continuation flag value for subsequent
c     input checking
   50 jflag=iflag
      kflag=0
c
c     save relerr and abserr for checking input on subsequent calls
      savre=relerr
      savae=abserr
c
c     restrict relative error tolerance to be at least as large as
c     2*eps+remin to avoid limiting precision difficulties arising
c     from impossible accuracy requests
c
      rer=twoeps+remin
      if (relerr .ge. rer) go to 55
c
c     relative error tolerance too small
      relerr=rer
      iflag=3
      kflag=3
      return
c
   55 dt=tout-t
c
      if (mflag .eq. 1) go to 60
      if (init .eq. 0) go to 65
      go to 80
c
c     initialization --
c                       set initialization completion indicator,init
c                       set indicator for too many output points,kop
c                       evaluate initial derivatives
c                       set counter for function evaluations,nfe
c                       evaluate initial derivatives
c                       set counter for function evaluations,nfe
c                       estimate starting stepsize
c
   60 init=0
      kop=0
c
      a=t
      call f(a,y,yp)
      nfe=1
      if (t .ne. tout) go to 65
      iflag=2
      return
c
c
   65 init=1
      h=dabs(dt)
      toln=0.
      do 70 k=1,neqn
        tol=relerr*dabs(y(k))+abserr
        if (tol .le. 0.) go to 70
        toln=tol
        ypk=dabs(yp(k))
        if (ypk*h**5 .gt. tol) h=(tol/ypk)**0.2d0
   70 continue
      if (toln .le. 0.0d0) h=0.0d0
      h=dmax1(h,u26*dmax1(dabs(t),dabs(dt)))
      jflag=isign(2,iflag)
c
c
c     set stepsize for integration in the direction from t to tout
c
   80 h=dsign(h,dt)
c
c     test to see if rkf45 is being severely impacted by too many
c     output points
c
      if (dabs(h) .ge. 2.0d0*dabs(dt)) kop=kop+1
      if (kop .ne. 100) go to 85
c
c     unnecessary frequency of output
      kop=0
      iflag=7
      return
c
   85 if (dabs(dt) .gt. u26*dabs(t)) go to 95
c
c     if too close to output point,extrapolate and return
c
      do 90 k=1,neqn
   90   y(k)=y(k)+dt*yp(k)
      a=tout
      call f(a,y,yp)
      nfe=nfe+1
      go to 300
c
c
c     initialize output point indicator
c
   95 output= .false.
c
c     to avoid premature underflow in the error tolerance function,
c     scale the error tolerances
c
      scale=2.0d0/relerr
      ae=scale*abserr
c
c
c     step by step integration
c
  100 hfaild= .false.
c
c     set smallest allowable stepsize
c
      hmin=u26*dabs(t)
c
c     adjust stepsize if necessary to hit the output point.
c     look ahead two steps to avoid drastic changes in the stepsize and
c     thus lessen the impact of output points on the code.
c
      dt=tout-t
      if (dabs(dt) .ge. 2.0d0*dabs(h)) go to 200
      if (dabs(dt) .gt. dabs(h)) go to 150
c
c     the next successful step will complete the integration to the
c     output point
c
      output= .true.
      h=dt
      go to 200
c
  150 h=0.5d0*dt
c
c
c
c     core integrator for taking a single step
c
c     the tolerances have been scaled to avoid premature underflow in
c     computing the error tolerance function et.
c     to avoid problems with zero crossings,relative error is measured
c     using the average of the magnitudes of the solution at the
c     beginning and end of a step.
c     the error estimate formula has been grouped to control loss of
c     significance.
c     to distinguish the various arguments, h is not permitted
c     to become smaller than 26 units of roundoff in t.
c     practical limits on the change in the stepsize are enforced to
c     smooth the stepsize selection process and to avoid excessive
c     chattering on problems having discontinuities.
c     to prevent unnecessary failures, the code uses 9/10 the stepsize
c     it estimates will succeed.
c     after a step failure, the stepsize is not allowed to increase for
c     the next attempted step. this makes the code more efficient on
c     problems having discontinuities and more effective in general
c     since local extrapolation is being used and extra caution seems
c     warranted.
c
c
c     test number of derivative function evaluations.
c     if okay,try to advance the integration from t to t+h
c
  200 if (nfe .le. maxnfe) go to 220
c
c     too much work
      iflag=4
      kflag=4
      return
c
c     advance an approximate solution over one step of length h
c
  220 call fehl(f,neqn,y,t,h,yp,f1,f2,f3,f4,f5,f1)
      nfe=nfe+5
c
c     compute and test allowable tolerances versus local error estimates
c     and remove scaling of tolerances. note that relative error is
c     measured with respect to the average of the magnitudes of the
c     solution at the beginning and end of the step.
c
      eeoet=0.0d0
      do 250 k=1,neqn
        et=dabs(y(k))+dabs(f1(k))+ae
        if (et .gt. 0.0d0) go to 240
c
c       inappropriate error tolerance
        iflag=5
        return
c
  240   ee=dabs((-2090.0d0*yp(k)+(21970.0d0*f3(k)-15048.0d0*f4(k)))+
     1                        (22528.0d0*f2(k)-27360.0d0*f5(k)))
  250   eeoet=dmax1(eeoet,ee/et)
c
      esttol=dabs(h)*eeoet*scale/752400.0d0
c
      if (esttol .le. 1.0d0) go to 260
c
c
c     unsuccessful step
c                       reduce the stepsize , try again
c                       the decrease is limited to a factor of 1/10
c
      hfaild= .true.
      output= .false.
      s=0.1d0
      if (esttol .lt. 59049.0d0) s=0.9d0/esttol**0.2d0
      h=s*h
      if (dabs(h) .gt. hmin) go to 200
c
c     requested error unattainable at smallest allowable stepsize
      iflag=6
      kflag=6
      return
c
c
c     successful step
c                        store solution at t+h
c                        and evaluate derivatives there
c
  260 t=t+h
      do 270 k=1,neqn
  270   y(k)=f1(k)
      a=t
      call f(a,y,yp)
      nfe=nfe+1
c
c
c                       choose next stepsize
c                       the increase is limited to a factor of 5
c                       if step failure has just occurred, next
c                          stepsize is not allowed to increase
c
      s=5.0d0
      if (esttol .gt. 1.889568d-4) s=0.9d0/esttol**0.2d0
      if (hfaild) s=dmin1(s,1.0d0)
      h=dsign(dmax1(s*dabs(h),hmin),h)
c
c     end of core integrator
c
c
c     should we take another step
c
      if (output) go to 300
      if (iflag .gt. 0) go to 100
c
c
c     integration successfully completed
c
c     one-step mode
      iflag=-2
      return
c
c     interval mode
  300 t=tout
      iflag=2
      return
c
      end
      subroutine fehl(f,neqn,y,t,h,yp,f1,f2,f3,f4,f5,s)
c
c     fehlberg fourth-fifth order runge-kutta method
c
c    fehl integrates a system of neqn first order
c    ordinary differential equations of the form
c             dy(i)/dt=f(t,y(1),---,y(neqn))
c    where the initial values y(i) and the initial derivatives
c    yp(i) are specified at the starting point t. fehl advances
c    the solution over the fixed step h and returns
c    the fifth order (sixth order accurate locally) solution
c    approximation at t+h in array s(i).
c    f1,---,f5 are arrays of dimension neqn which are needed
c    for internal storage.
c    the formulas have been grouped to control loss of significance.
c    fehl should be called with an h not smaller than 13 units of
c    roundoff in t so that the various independent arguments can be
c    distinguished.
c
c
      integer  neqn
      double precision  y(neqn),t,h,yp(neqn),f1(neqn),f2(neqn),
     1  f3(neqn),f4(neqn),f5(neqn),s(neqn)
c
      double precision  ch
      integer  k
c
      ch=h/4.0d0
      do 221 k=1,neqn
  221   f5(k)=y(k)+ch*yp(k)
      call f(t+ch,f5,f1)
c
      ch=3.0d0*h/32.0d0
      do 222 k=1,neqn
  222   f5(k)=y(k)+ch*(yp(k)+3.0d0*f1(k))
      call f(t+3.0d0*h/8.0d0,f5,f2)
c
      ch=h/2197.0d0
      do 223 k=1,neqn
  223   f5(k)=y(k)+ch*(1932.0d0*yp(k)+(7296.0d0*f2(k)-7200.0d0*f1(k)))
      call f(t+12.0d0*h/13.0d0,f5,f3)
c
      ch=h/4104.0d0
      do 224 k=1,neqn
  224   f5(k)=y(k)+ch*((8341.0d0*yp(k)-845.0d0*f3(k))+
     1                            (29440.0d0*f2(k)-32832.0d0*f1(k)))
      call f(t+h,f5,f4)
c
      ch=h/20520.0d0
      do 225 k=1,neqn
  225   f1(k)=y(k)+ch*((-6080.0d0*yp(k)+(9295.0d0*f3(k)-
     1         5643.0d0*f4(k)))+(41040.0d0*f1(k)-28352.0d0*f2(k)))
      call f(t+h/2.0d0,f1,f5)
c
c     compute approximate solution at t+h
c
      ch=h/7618050.0d0
      do 230 k=1,neqn
  230   s(k)=y(k)+ch*((902880.0d0*yp(k)+(3855735.0d0*f3(k)-
     1        1371249.0d0*f4(k)))+(3953664.0d0*f2(k)+
     2        277020.0d0*f5(k)))
c
      return
      end
      function d1mach ( i )

c*********************************************************************72
c
cc D1MACH returns double precision machine-dependent constants.
c
c  Discussion:
c
c    D1MACH can be used to obtain machine-dependent parameters
c    for the local machine environment.  It is a function
c    with one input argument, and can be called as follows:
c
c      D = D1MACH ( I )
c
c    where I=1,...,5.  The output value of D above is
c    determined by the input value of I:.
c
c    D1MACH ( 1) = B^(EMIN-1), the smallest positive magnitude.
c    D1MACH ( 2) = B^EMAX*(1 - B^(-T)), the largest magnitude.
c    D1MACH ( 3) = B^(-T), the smallest relative spacing.
c    D1MACH ( 4) = B^(1-T), the largest relative spacing.
c    D1MACH ( 5) = LOG10(B)
c
c  Modified:
c
c    06 December 2006
c
c  Author:
c
c    Phyllis Fox, Andrew Hall, Norman Schryer
c
c  Reference:
c
c    Phyllis Fox, Andrew Hall, Norman Schryer,
c    Algorithm 528:
c    Framework for a Portable Library,
c    ACM Transactions on Mathematical Software,
c    Volume 4, Number 2, June 1978, page 176-188.
c
c  Parameters:
c
c    Input, integer I, the index of the desired constant.
c
c    Output, double precision D1MACH, the value of the constant.
c
      implicit none

      double precision d1mach
      integer diver(4)
      double precision dmach(5)
      integer i
      integer large(4)
      integer log10(4)
      integer right(4)
      integer small(4)

      equivalence ( dmach(1), small(1) )
      equivalence ( dmach(2), large(1) )
      equivalence ( dmach(3), right(1) )
      equivalence ( dmach(4), diver(1) )
      equivalence ( dmach(5), log10(1) )
c
c     MACHINE CONSTANTS FOR IEEE ARITHMETIC MACHINES, SUCH AS THE AT&T
c     3B SERIES AND MOTOROLA 68000 BASED MACHINES (E.G. SUN 3 AND AT&T
c     PC 7300), IN WHICH THE MOST SIGNIFICANT BYTE IS STORED FIRST.
c
c === MACHINE = IEEE.MOST-SIG-BYTE-FIRST
c === MACHINE = SUN
c === MACHINE = 68000
c === MACHINE = ATT.3B
c === MACHINE = ATT.7300
c     DATA SMALL(1),SMALL(2) /    1048576,          0 /
c     DATA LARGE(1),LARGE(2) / 2146435071,         -1 /
c     DATA RIGHT(1),RIGHT(2) / 1017118720,          0 /
c     DATA DIVER(1),DIVER(2) / 1018167296,          0 /
c     DATA LOG10(1),LOG10(2) / 1070810131, 1352628735 /
c
c     MACHINE CONSTANTS FOR IEEE ARITHMETIC MACHINES AND 8087-BASED
c     MICROS, SUCH AS THE IBM PC AND AT&T 6300, IN WHICH THE LEAST
c     SIGNIFICANT BYTE IS STORED FIRST.
c
c === MACHINE = IEEE.LEAST-SIG-BYTE-FIRST
c === MACHINE = 8087
c === MACHINE = IBM.PC
c === MACHINE = ATT.6300
c
       data small(1),small(2) /          0,    1048576 /
       data large(1),large(2) /         -1, 2146435071 /
       data right(1),right(2) /          0, 1017118720 /
       data diver(1),diver(2) /          0, 1018167296 /
       data log10(1),log10(2) / 1352628735, 1070810131 /
c
c     MACHINE CONSTANTS FOR AMDAHL MACHINES.
c
c === MACHINE = AMDAHL
c      DATA SMALL(1),SMALL(2) /    1048576,          0 /
c      DATA LARGE(1),LARGE(2) / 2147483647,         -1 /
c      DATA RIGHT(1),RIGHT(2) /  856686592,          0 /
c      DATA DIVER(1),DIVER(2) /  873463808,          0 /
c      DATA LOG10(1),LOG10(2) / 1091781651, 1352628735 /
c
c     MACHINE CONSTANTS FOR THE BURROUGHS 1700 SYSTEM.
c
c === MACHINE = BURROUGHS.1700
c      DATA SMALL(1) / ZC00800000 /
c      DATA SMALL(2) / Z000000000 /
c      DATA LARGE(1) / ZDFFFFFFFF /
c      DATA LARGE(2) / ZFFFFFFFFF /
c      DATA RIGHT(1) / ZCC5800000 /
c      DATA RIGHT(2) / Z000000000 /
c      DATA DIVER(1) / ZCC6800000 /
c      DATA DIVER(2) / Z000000000 /
c      DATA LOG10(1) / ZD00E730E7 /
c      DATA LOG10(2) / ZC77800DC0 /
c
c     MACHINE CONSTANTS FOR THE BURROUGHS 5700 SYSTEM.
c
c === MACHINE = BURROUGHS.5700
c      DATA SMALL(1) / O1771000000000000 /
c      DATA SMALL(2) / O0000000000000000 /
c      DATA LARGE(1) / O0777777777777777 /
c      DATA LARGE(2) / O0007777777777777 /
c      DATA RIGHT(1) / O1461000000000000 /
c      DATA RIGHT(2) / O0000000000000000 /
c      DATA DIVER(1) / O1451000000000000 /
c      DATA DIVER(2) / O0000000000000000 /
c      DATA LOG10(1) / O1157163034761674 /
c      DATA LOG10(2) / O0006677466732724 /
c
c     MACHINE CONSTANTS FOR THE BURROUGHS 6700/7700 SYSTEMS.
c
c === MACHINE = BURROUGHS.6700
c === MACHINE = BURROUGHS.7700
c      DATA SMALL(1) / O1771000000000000 /
c      DATA SMALL(2) / O7770000000000000 /
c      DATA LARGE(1) / O0777777777777777 /
c      DATA LARGE(2) / O7777777777777777 /
c      DATA RIGHT(1) / O1461000000000000 /
c      DATA RIGHT(2) / O0000000000000000 /
c      DATA DIVER(1) / O1451000000000000 /
c      DATA DIVER(2) / O0000000000000000 /
c      DATA LOG10(1) / O1157163034761674 /
c      DATA LOG10(2) / O0006677466732724 /
c
c     MACHINE CONSTANTS FOR THE CONVEX C-120 (NATIVE MODE)
c     WITH OR WITHOUT -R8 OPTION
c
c === MACHINE = CONVEX.C1
c === MACHINE = CONVEX.C1.R8
c      DATA DMACH(1) / 5.562684646268007D-309 /
c      DATA DMACH(2) / 8.988465674311577D+307 /
c      DATA DMACH(3) / 1.110223024625157D-016 /
c      DATA DMACH(4) / 2.220446049250313D-016 /
c      DATA DMACH(5) / 3.010299956639812D-001 /
c
c     MACHINE CONSTANTS FOR THE CONVEX C-120 (IEEE MODE)
c     WITH OR WITHOUT -R8 OPTION
c
c === MACHINE = CONVEX.C1.IEEE
c === MACHINE = CONVEX.C1.IEEE.R8
c      DATA DMACH(1) / 2.225073858507202D-308 /
c      DATA DMACH(2) / 1.797693134862315D+308 /
c      DATA DMACH(3) / 1.110223024625157D-016 /
c      DATA DMACH(4) / 2.220446049250313D-016 /
c      DATA DMACH(5) / 3.010299956639812D-001 /
c
c     MACHINE CONSTANTS FOR THE CYBER 170/180 SERIES USING NOS (FTN5).
c
c === MACHINE = CYBER.170.NOS
c === MACHINE = CYBER.180.NOS
c      DATA SMALL(1) / O"00604000000000000000" /
c      DATA SMALL(2) / O"00000000000000000000" /
c      DATA LARGE(1) / O"37767777777777777777" /
c      DATA LARGE(2) / O"37167777777777777777" /
c      DATA RIGHT(1) / O"15604000000000000000" /
c      DATA RIGHT(2) / O"15000000000000000000" /
c      DATA DIVER(1) / O"15614000000000000000" /
c      DATA DIVER(2) / O"15010000000000000000" /
c      DATA LOG10(1) / O"17164642023241175717" /
c      DATA LOG10(2) / O"16367571421742254654" /
c
c     MACHINE CONSTANTS FOR THE CDC 180 SERIES USING NOS/VE
c
c === MACHINE = CYBER.180.NOS/VE
c      DATA SMALL(1) / Z"3001800000000000" /
c      DATA SMALL(2) / Z"3001000000000000" /
c      DATA LARGE(1) / Z"4FFEFFFFFFFFFFFE" /
c      DATA LARGE(2) / Z"4FFE000000000000" /
c      DATA RIGHT(1) / Z"3FD2800000000000" /
c      DATA RIGHT(2) / Z"3FD2000000000000" /
c      DATA DIVER(1) / Z"3FD3800000000000" /
c      DATA DIVER(2) / Z"3FD3000000000000" /
c      DATA LOG10(1) / Z"3FFF9A209A84FBCF" /
c      DATA LOG10(2) / Z"3FFFF7988F8959AC" /
c
c     MACHINE CONSTANTS FOR THE CYBER 205
c
c === MACHINE = CYBER.205
c      DATA SMALL(1) / X'9000400000000000' /
c      DATA SMALL(2) / X'8FD1000000000000' /
c      DATA LARGE(1) / X'6FFF7FFFFFFFFFFF' /
c      DATA LARGE(2) / X'6FD07FFFFFFFFFFF' /
c      DATA RIGHT(1) / X'FF74400000000000' /
c      DATA RIGHT(2) / X'FF45000000000000' /
c      DATA DIVER(1) / X'FF75400000000000' /
c      DATA DIVER(2) / X'FF46000000000000' /
c      DATA LOG10(1) / X'FFD04D104D427DE7' /
c      DATA LOG10(2) / X'FFA17DE623E2566A' /
c
c     MACHINE CONSTANTS FOR THE CDC 6000/7000 SERIES.
c
c === MACHINE = CDC.6000
c === MACHINE = CDC.7000
c      DATA SMALL(1) / 00604000000000000000B /
c      DATA SMALL(2) / 00000000000000000000B /
c      DATA LARGE(1) / 37767777777777777777B /
c      DATA LARGE(2) / 37167777777777777777B /
c      DATA RIGHT(1) / 15604000000000000000B /
c      DATA RIGHT(2) / 15000000000000000000B /
c      DATA DIVER(1) / 15614000000000000000B /
c      DATA DIVER(2) / 15010000000000000000B /
c      DATA LOG10(1) / 17164642023241175717B /
c      DATA LOG10(2) / 16367571421742254654B /
c
c     MACHINE CONSTANTS FOR THE CRAY 1, XMP, 2, AND 3.
c
c === MACHINE = CRAY
c      DATA SMALL(1) / 201354000000000000000B /
c      DATA SMALL(2) / 000000000000000000000B /
c      DATA LARGE(1) / 577767777777777777777B /
c      DATA LARGE(2) / 000007777777777777776B /
c      DATA RIGHT(1) / 376434000000000000000B /
c      DATA RIGHT(2) / 000000000000000000000B /
c      DATA DIVER(1) / 376444000000000000000B /
c      DATA DIVER(2) / 000000000000000000000B /
c      DATA LOG10(1) / 377774642023241175717B /
c      DATA LOG10(2) / 000007571421742254654B /
c
c     MACHINE CONSTANTS FOR THE DATA GENERAL ECLIPSE S/200
c
c     NOTE - IT MAY BE APPROPRIATE TO INCLUDE THE FOLLOWING LINE -
c     STATIC DMACH(5)
c
c === MACHINE = DATA_GENERAL.ECLIPSE.S/200
c      DATA SMALL/20K,3*0/,LARGE/77777K,3*177777K/
c      DATA RIGHT/31420K,3*0/,DIVER/32020K,3*0/
c      DATA LOG10/40423K,42023K,50237K,74776K/
c
c     ELXSI 6400
c
c === MACHINE = ELSXI.6400
c      DATA SMALL(1), SMALL(2) / '00100000'X,'00000000'X /
c      DATA LARGE(1), LARGE(2) / '7FEFFFFF'X,'FFFFFFFF'X /
c      DATA RIGHT(1), RIGHT(2) / '3CB00000'X,'00000000'X /
c      DATA DIVER(1), DIVER(2) / '3CC00000'X,'00000000'X /
c      DATA LOG10(1), DIVER(2) / '3FD34413'X,'509F79FF'X /
c
c     MACHINE CONSTANTS FOR THE HARRIS 220
c     MACHINE CONSTANTS FOR THE HARRIS SLASH 6 AND SLASH 7
c
c === MACHINE = HARRIS.220
c === MACHINE = HARRIS.SLASH6
c === MACHINE = HARRIS.SLASH7
c      DATA SMALL(1),SMALL(2) / '20000000, '00000201 /
c      DATA LARGE(1),LARGE(2) / '37777777, '37777577 /
c      DATA RIGHT(1),RIGHT(2) / '20000000, '00000333 /
c      DATA DIVER(1),DIVER(2) / '20000000, '00000334 /
c      DATA LOG10(1),LOG10(2) / '23210115, '10237777 /
c
c     MACHINE CONSTANTS FOR THE HONEYWELL 600/6000 SERIES.
c     MACHINE CONSTANTS FOR THE HONEYWELL DPS 8/70 SERIES.
c
c === MACHINE = HONEYWELL.600/6000
c === MACHINE = HONEYWELL.DPS.8/70
c      DATA SMALL(1),SMALL(2) / O402400000000, O000000000000 /
c      DATA LARGE(1),LARGE(2) / O376777777777, O777777777777 /
c      DATA RIGHT(1),RIGHT(2) / O604400000000, O000000000000 /
c      DATA DIVER(1),DIVER(2) / O606400000000, O000000000000 /
c      DATA LOG10(1),LOG10(2) / O776464202324, O117571775714 /
c
c      MACHINE CONSTANTS FOR THE HP 2100
c      3 WORD DOUBLE PRECISION OPTION WITH FTN4
c
c === MACHINE = HP.2100.3_WORD_DP
c      DATA SMALL(1), SMALL(2), SMALL(3) / 40000B,       0,       1 /
c      DATA LARGE(1), LARGE(2), LARGE(3) / 77777B, 177777B, 177776B /
c      DATA RIGHT(1), RIGHT(2), RIGHT(3) / 40000B,       0,    265B /
c      DATA DIVER(1), DIVER(2), DIVER(3) / 40000B,       0,    276B /
c      DATA LOG10(1), LOG10(2), LOG10(3) / 46420B,  46502B,  77777B /
c
c      MACHINE CONSTANTS FOR THE HP 2100
c      4 WORD DOUBLE PRECISION OPTION WITH FTN4
c
c === MACHINE = HP.2100.4_WORD_DP
c      DATA SMALL(1), SMALL(2) /  40000B,       0 /
c      DATA SMALL(3), SMALL(4) /       0,       1 /
c      DATA LARGE(1), LARGE(2) /  77777B, 177777B /
c      DATA LARGE(3), LARGE(4) / 177777B, 177776B /
c      DATA RIGHT(1), RIGHT(2) /  40000B,       0 /
c      DATA RIGHT(3), RIGHT(4) /       0,    225B /
c      DATA DIVER(1), DIVER(2) /  40000B,       0 /
c      DATA DIVER(3), DIVER(4) /       0,    227B /
c      DATA LOG10(1), LOG10(2) /  46420B,  46502B /
c      DATA LOG10(3), LOG10(4) /  76747B, 176377B /
c
c     HP 9000
c
c      D1MACH(1) = 2.8480954D-306
c      D1MACH(2) = 1.40444776D+306
c      D1MACH(3) = 2.22044605D-16
c      D1MACH(4) = 4.44089210D-16
c      D1MACH(5) = 3.01029996D-1
c
c === MACHINE = HP.9000
c      DATA SMALL(1), SMALL(2) / 00040000000B, 00000000000B /
c      DATA LARGE(1), LARGE(2) / 17737777777B, 37777777777B /
c      DATA RIGHT(1), RIGHT(2) / 07454000000B, 00000000000B /
c      DATA DIVER(1), DIVER(2) / 07460000000B, 00000000000B /
c      DATA LOG10(1), LOG10(2) / 07764642023B, 12047674777B /
c
c     MACHINE CONSTANTS FOR THE IBM 360/370 SERIES,
c     THE XEROX SIGMA 5/7/9, THE SEL SYSTEMS 85/86, AND
c     THE INTERDATA 3230 AND INTERDATA 7/32.
c
c === MACHINE = IBM.360
c === MACHINE = IBM.370
c === MACHINE = XEROX.SIGMA.5
c === MACHINE = XEROX.SIGMA.7
c === MACHINE = XEROX.SIGMA.9
c === MACHINE = SEL.85
c === MACHINE = SEL.86
c === MACHINE = INTERDATA.3230
c === MACHINE = INTERDATA.7/32
c      DATA SMALL(1),SMALL(2) / Z00100000, Z00000000 /
c      DATA LARGE(1),LARGE(2) / Z7FFFFFFF, ZFFFFFFFF /
c      DATA RIGHT(1),RIGHT(2) / Z33100000, Z00000000 /
c      DATA DIVER(1),DIVER(2) / Z34100000, Z00000000 /
c      DATA LOG10(1),LOG10(2) / Z41134413, Z509F79FF /
c
c     MACHINE CONSTANTS FOR THE INTERDATA 8/32
c     WITH THE UNIX SYSTEM FORTRAN 77 COMPILER.
c
c     FOR THE INTERDATA FORTRAN VII COMPILER REPLACE
c     THE Z'S SPECIFYING HEX CONSTANTS WITH Y'S.
c
c === MACHINE = INTERDATA.8/32.UNIX
c      DATA SMALL(1),SMALL(2) / Z'00100000', Z'00000000' /
c      DATA LARGE(1),LARGE(2) / Z'7EFFFFFF', Z'FFFFFFFF' /
c      DATA RIGHT(1),RIGHT(2) / Z'33100000', Z'00000000' /
c      DATA DIVER(1),DIVER(2) / Z'34100000', Z'00000000' /
c      DATA LOG10(1),LOG10(2) / Z'41134413', Z'509F79FF' /
c
c     MACHINE CONSTANTS FOR THE PDP-10 (KA PROCESSOR).
c
c === MACHINE = PDP-10.KA
c      DATA SMALL(1),SMALL(2) / "033400000000, "000000000000 /
c      DATA LARGE(1),LARGE(2) / "377777777777, "344777777777 /
c      DATA RIGHT(1),RIGHT(2) / "113400000000, "000000000000 /
c      DATA DIVER(1),DIVER(2) / "114400000000, "000000000000 /
c      DATA LOG10(1),LOG10(2) / "177464202324, "144117571776 /
c
c     MACHINE CONSTANTS FOR THE PDP-10 (KI PROCESSOR).
c
c === MACHINE = PDP-10.KI
c      DATA SMALL(1),SMALL(2) / "000400000000, "000000000000 /
c      DATA LARGE(1),LARGE(2) / "377777777777, "377777777777 /
c      DATA RIGHT(1),RIGHT(2) / "103400000000, "000000000000 /
c      DATA DIVER(1),DIVER(2) / "104400000000, "000000000000 /
c      DATA LOG10(1),LOG10(2) / "177464202324, "047674776746 /
c
c     MACHINE CONSTANTS FOR PDP-11 FORTRAN SUPPORTING
c     32-BIT INTEGERS (EXPRESSED IN INTEGER AND OCTAL).
c
c === MACHINE = PDP-11.32-BIT
c      DATA SMALL(1),SMALL(2) /    8388608,           0 /
c      DATA LARGE(1),LARGE(2) / 2147483647,          -1 /
c      DATA RIGHT(1),RIGHT(2) /  612368384,           0 /
c      DATA DIVER(1),DIVER(2) /  620756992,           0 /
c      DATA LOG10(1),LOG10(2) / 1067065498, -2063872008 /
c
c      DATA SMALL(1),SMALL(2) / O00040000000, O00000000000 /
c      DATA LARGE(1),LARGE(2) / O17777777777, O37777777777 /
c      DATA RIGHT(1),RIGHT(2) / O04440000000, O00000000000 /
c      DATA DIVER(1),DIVER(2) / O04500000000, O00000000000 /
c      DATA LOG10(1),LOG10(2) / O07746420232, O20476747770 /
c
c     MACHINE CONSTANTS FOR PDP-11 FORTRAN SUPPORTING
c     16-BIT INTEGERS (EXPRESSED IN INTEGER AND OCTAL).
c
c === MACHINE = PDP-11.16-BIT
c      DATA SMALL(1),SMALL(2) /    128,      0 /
c      DATA SMALL(3),SMALL(4) /      0,      0 /
c      DATA LARGE(1),LARGE(2) /  32767,     -1 /
c      DATA LARGE(3),LARGE(4) /     -1,     -1 /
c      DATA RIGHT(1),RIGHT(2) /   9344,      0 /
c      DATA RIGHT(3),RIGHT(4) /      0,      0 /
c      DATA DIVER(1),DIVER(2) /   9472,      0 /
c      DATA DIVER(3),DIVER(4) /      0,      0 /
c      DATA LOG10(1),LOG10(2) /  16282,   8346 /
c      DATA LOG10(3),LOG10(4) / -31493, -12296 /
c
c      DATA SMALL(1),SMALL(2) / O000200, O000000 /
c      DATA SMALL(3),SMALL(4) / O000000, O000000 /
c      DATA LARGE(1),LARGE(2) / O077777, O177777 /
c      DATA LARGE(3),LARGE(4) / O177777, O177777 /
c      DATA RIGHT(1),RIGHT(2) / O022200, O000000 /
c      DATA RIGHT(3),RIGHT(4) / O000000, O000000 /
c      DATA DIVER(1),DIVER(2) / O022400, O000000 /
c      DATA DIVER(3),DIVER(4) / O000000, O000000 /
c      DATA LOG10(1),LOG10(2) / O037632, O020232 /
c      DATA LOG10(3),LOG10(4) / O102373, O147770 /
c
c     MACHINE CONSTANTS FOR THE SEQUENT BALANCE 8000
c
c === MACHINE = SEQUENT.BALANCE.8000
c      DATA SMALL(1),SMALL(2) / $00000000,  $00100000 /
c      DATA LARGE(1),LARGE(2) / $FFFFFFFF,  $7FEFFFFF /
c      DATA RIGHT(1),RIGHT(2) / $00000000,  $3CA00000 /
c      DATA DIVER(1),DIVER(2) / $00000000,  $3CB00000 /
c      DATA LOG10(1),LOG10(2) / $509F79FF,  $3FD34413 /
c
c     MACHINE CONSTANTS FOR THE UNIVAC 1100 SERIES. FTN COMPILER
c
c === MACHINE = UNIVAC.1100
c      DATA SMALL(1),SMALL(2) / O000040000000, O000000000000 /
c      DATA LARGE(1),LARGE(2) / O377777777777, O777777777777 /
c      DATA RIGHT(1),RIGHT(2) / O170540000000, O000000000000 /
c      DATA DIVER(1),DIVER(2) / O170640000000, O000000000000 /
c      DATA LOG10(1),LOG10(2) / O177746420232, O411757177572 /
c
c     MACHINE CONSTANTS FOR VAX 11/780
c     (EXPRESSED IN INTEGER AND HEXADECIMAL)
c    *** THE INTEGER FORMAT SHOULD BE OK FOR UNIX SYSTEMS***
c
c === MACHINE = VAX.11/780
c      DATA SMALL(1), SMALL(2) /        128,           0 /
c      DATA LARGE(1), LARGE(2) /     -32769,          -1 /
c      DATA RIGHT(1), RIGHT(2) /       9344,           0 /
c      DATA DIVER(1), DIVER(2) /       9472,           0 /
c      DATA LOG10(1), LOG10(2) /  546979738,  -805796613 /
c
c    ***THE HEX FORMAT BELOW MAY NOT BE SUITABLE FOR UNIX SYSYEMS***
c      DATA SMALL(1), SMALL(2) / Z00000080, Z00000000 /
c      DATA LARGE(1), LARGE(2) / ZFFFF7FFF, ZFFFFFFFF /
c      DATA RIGHT(1), RIGHT(2) / Z00002480, Z00000000 /
c      DATA DIVER(1), DIVER(2) / Z00002500, Z00000000 /
c      DATA LOG10(1), LOG10(2) / Z209A3F9A, ZCFF884FB /
c
c   MACHINE CONSTANTS FOR VAX 11/780 (G-FLOATING)
c     (EXPRESSED IN INTEGER AND HEXADECIMAL)
c    *** THE INTEGER FORMAT SHOULD BE OK FOR UNIX SYSTEMS***
c
c      DATA SMALL(1), SMALL(2) /         16,           0 /
c      DATA LARGE(1), LARGE(2) /     -32769,          -1 /
c      DATA RIGHT(1), RIGHT(2) /      15552,           0 /
c      DATA DIVER(1), DIVER(2) /      15568,           0 /
c      DATA LOG10(1), LOG10(2) /  1142112243, 2046775455 /
c
c    ***THE HEX FORMAT BELOW MAY NOT BE SUITABLE FOR UNIX SYSYEMS***
c      DATA SMALL(1), SMALL(2) / Z00000010, Z00000000 /
c      DATA LARGE(1), LARGE(2) / ZFFFF7FFF, ZFFFFFFFF /
c      DATA RIGHT(1), RIGHT(2) / Z00003CC0, Z00000000 /
c      DATA DIVER(1), DIVER(2) / Z00003CD0, Z00000000 /
c      DATA LOG10(1), LOG10(2) / Z44133FF3, Z79FF509F /
c
      if ( i .lt. 1  .or.  5 .lt. i ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'D1MACH - Fatal error!'
        write ( *, '(a)' ) '  I out of bounds.'
        stop
      end if

      d1mach = dmach(i)

      return
      end

