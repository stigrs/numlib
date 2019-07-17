
/* The MIT License (MIT)
 *
 * Copyright (c) 2014 ESSS
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

/* Modifications by Stig Rune Sellevag, 22 Jan 2019 */

#ifndef NUMLIB_MATH_CQUADPACK_H
#define NUMLIB_MATH_CQUADPACK_H

#define _USE_MATH_DEFINES

#include <math.h>
#include <float.h>

#define uflow_ DBL_MIN
#define oflow_ DBL_MAX
#define epmach DBL_EPSILON
#define LIMIT 500

#ifndef FALSE
#define FALSE 0
#endif
#ifndef TRUE
#define TRUE 1
#endif

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

typedef double(*dq_function_type)(double, void*);

/* Integration routines */
/* Gauss-Kronrod for integration over finite range. */
double cquadpack_G_K21(dq_function_type f,
                       double a,
                       double b,
                       double* abserr,
                       double* resabs,
                       double* resasc,
                       void* user_data);

/* Gauss-Kronrod for integration over infinite range. */
double cquadpack_G_K15I(dq_function_type f,
                        double boun,
                        int inf,
                        double a,
                        double b,
                        double* abserr,
                        double* resabs,
                        double* resasc,
                        void* user_data);

double cquadpack_dqext(
    int* n, double epstab[], double* abserr, double res3la[], int* nres);

void cquadpack_dqsort(int limit,
                      int last,
                      int* maxerr,
                      double* ermax,
                      double elist[],
                      int iord[],
                      int* nrmax);

double cquadpack_dqagi(dq_function_type f,
                       double bound,
                       int inf,
                       double epsabs,
                       double epsrel,
                       double* abserr,
                       int* neval,
                       int* ier,
                       void* user_data);

double cquadpack_dqags(dq_function_type f,
                       double a,
                       double b,
                       double epsabs,
                       double epsrel,
                       double* abserr,
                       int* neval,
                       int* ier,
                       void* user_data);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* NUMLIB_MATH_CQUADPACK_H */
