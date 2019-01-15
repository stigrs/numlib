/***
 *        Created:  2018-08-14
 *
 *         Author:  Dilawar Singh <dilawars@ncbs.res.in>
 *   Organization:  NCBS Bangalore
 *        License:  MIT License
 */

#ifndef NUMLIB_MATH_LSODA_H
#define NUMLIB_MATH_LSODA_H

#include <array>
#include <cmath>
#include <memory>
#include <vector>

/* --------------------------------------------------------------------------*/
/**
 * @Synopsis  Type definition of LSODA ode system. See the file test_LSODA.cpp
 * for an example.
 *
 * @Param time, double
 * @Param y, array of double.
 * @Param dydt, array of double
 * @Param data, void*
 *
 * @Returns void
 */
/* ----------------------------------------------------------------------------*/
typedef void (*Lsoda_sys_t)(double t, double* y, double* dydt, void*);

namespace Numlib {

class Lsoda {

public:
    Lsoda() = delete;

    Lsoda(Lsoda_sys_t f, double atol = 1.0e-6, double rtol = 1.0e-6);

    ~Lsoda() = default;

    void
    integrate(Vec<double>& y, double& t, double tout, void* data = nullptr);

private:
    std::size_t idamax1(const std::vector<double>& dx,
                        const std::size_t n,
                        const std::size_t offset);

    void dscal1(const double da,
                std::vector<double>& dx,
                const std::size_t n,
                const std::size_t offset);

    double ddot1(const std::vector<double>& a,
                 const std::vector<double>& b,
                 const std::size_t n,
                 const std::size_t offsetA,
                 const std::size_t offsetB);

    void daxpy1(const double da,
                const std::vector<double>& dx,
                std::vector<double>& dy,
                const std::size_t n,
                const std::size_t offsetX,
                const std::size_t offsetY);

    void dgesl(const std::vector<std::vector<double>>& a,
               const std::size_t n,
               std::vector<int>& ipvt,
               std::vector<double>& b,
               const std::size_t job);

    void dgefa(std::vector<std::vector<double>>& a,
               const std::size_t n,
               std::vector<int>& ipvt,
               std::size_t* const info);

    void prja(const std::size_t neq,
              std::vector<double>& y,
              Lsoda_sys_t f,
              void* _data);

    void lsoda(Lsoda_sys_t f,
               const std::size_t neq,
               std::vector<double>& y,
               double* t,
               double tout,
               int itask,
               int* istate,
               int iopt,
               int jt,
               std::array<int, 7>& iworks,
               std::array<double, 4>& rworks,
               void* _data);

    void correction(const std::size_t neq,
                    std::vector<double>& y,
                    Lsoda_sys_t f,
                    std::size_t* corflag,
                    double pnorm,
                    double* del,
                    double* delp,
                    double* told,
                    std::size_t* ncf,
                    double* rh,
                    std::size_t* m,
                    void* _data);

    void stoda(const std::size_t neq,
               std::vector<double>& y,
               Lsoda_sys_t f,
               void* _data);

    void terminate(int* istate);
    void terminate2(std::vector<double>& y, double* t);
    void successreturn(std::vector<double>& y,
                       double* t,
                       int itask,
                       int ihit,
                       double tcrit,
                       int* istate);
    void ewset(const std::vector<double>& ycur);
    void resetcoeff(void);
    void solsy(std::vector<double>& y);
    void endstoda(void);
    void orderswitch(double* rhup,
                     double dsm,
                     double* pdh,
                     double* rh,
                     std::size_t* orderflag);
    void intdy(double t, int k, std::vector<double>& dky, int* iflag);
    void corfailure(double* told,
                    double* rh,
                    std::size_t* ncf,
                    std::size_t* corflag);
    void methodswitch(double dsm, double pnorm, double* pdh, double* rh);
    void cfode(int meth_);
    void scaleh(double* rh, double* pdh);
    double fnorm(int n,
                 const std::vector<std::vector<double>>& a,
                 const std::vector<double>& w);
    double vmnorm(const std::size_t n,
                  const std::vector<double>& v,
                  const std::vector<double>& w);

    static bool abs_compare(double a, double b);

    Lsoda_sys_t fsys;

    std::size_t ml;
    std::size_t mu;
    std::size_t imxer;
    double sqrteta;

    // NOTE: initialize in default constructor. Older compiler e.g. 4.8.4 would
    // produce error if these are initialized here. With newer compiler,
    // initialization can be done here.
    std::array<std::size_t, 3> mord;
    std::array<double, 13> sm1;

    std::array<double, 14> el;  // = {0};
    std::array<double, 13> cm1; // = {0};
    std::array<double, 6> cm2;  // = {0};

    std::array<std::array<double, 14>, 13> elco;
    std::array<std::array<double, 4>, 13> tesco;

    std::size_t illin;
    std::size_t init;
    std::size_t ierpj;
    std::size_t iersl;
    std::size_t jcur;
    std::size_t l;
    std::size_t miter;
    std::size_t maxord;
    std::size_t maxcor;
    std::size_t msbp;
    std::size_t mxncf;

    int kflag;
    int jstart;

    std::size_t ixpr = 0;
    std::size_t jtyp;
    std::size_t mused;
    std::size_t mxordn;
    std::size_t mxords = 12;
    std::size_t meth_;

    std::size_t n;
    std::size_t nq;
    std::size_t nst;
    std::size_t nfe;
    std::size_t nje;
    std::size_t nqu;
    std::size_t mxstep;
    std::size_t mxhnil;
    std::size_t nslast;
    std::size_t nhnil;
    std::size_t ntrep;
    std::size_t nyh;

    double ccmax;
    double el0;
    double h_ = .0;
    double hmin;
    double hmxi;
    double hu;
    double rc;
    double tn_ = 0.0;
    double tsw;
    double pdnorm;
    double conit;
    double crate;
    double hold;
    double rmax;

    std::size_t ialth;
    std::size_t ipup;
    std::size_t lmax;
    std::size_t nslp;
    double pdest;
    double pdlast;
    double ratio;
    int icount;
    int irflag;

    std::vector<double> ewt;
    std::vector<double> savf;
    std::vector<double> acor;
    std::vector<std::vector<double>> yh_;
    std::vector<std::vector<double>> wm_;

    std::vector<int> ipvt;

    int itol_ = 2;
    std::vector<double> rtol_;
    std::vector<double> atol_;
};

} // namespace Numlib

#endif // NUMLIB_MATH_LSODA_H

