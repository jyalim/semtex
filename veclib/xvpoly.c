/*****************************************************************************
 * xvpoly: Hoerner polynomial evaluation of a vector.
 *
 * m is the order of the polynomial and its coefficients are stored in c in
 * descending order: c[0] is the constant term.
 *
 * $Id: xvpoly.c,v 8.1 2015/04/20 11:14:21 hmb Exp $
 *****************************************************************************/

#include <cfemdef.h>

#if defined(__uxp__)
#pragma global novrec
#pragma global noalias
#endif

void dvpoly (int_t n,
	     const double* x, int_t incx, int_t m,
	     const double* c, int_t incc, 
	           double* y, int_t incy)
{
  register int_t i, j;
  register double  sum, xval;
  const    double  *csave, *cp;

  csave = c;
  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incx : 0;

  for (i = 0; i < n; i++) {
    c    = (incc<0) ? csave+(-n+1)*incc : csave;
    cp   = c + incc;
    sum  = c[0];
    xval = x[i*incx];
    for (j = 0; j < m; j++) sum = sum * xval + cp[j*incc];
    y[i*incy] = sum;
  }
}


void svpoly (int_t n,
	     const float* x, int_t incx, int_t m,
	     const float* c, int_t incc, 
	           float* y, int_t incy)
{
  register int_t i, j;
  register float   sum, xval;
  const    float   *csave, *cp;

  csave = c;
  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incx : 0;

  for (i = 0; i < n; i++) {
    c    = (incc<0) ? csave+(-n+1)*incc : csave;
    cp   = c + incc;
    sum  = c[0];
    xval = x[i*incx];
    for (j = 0; j < m; j++) sum = sum * xval + cp[j*incc];
    y[i*incy] = sum;
  }
}
