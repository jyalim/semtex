/*****************************************************************************
 * xvamax:  z[i] = MAX(ABS(x[i]), ABS(y[i])).
 *
 * $Id: xvamax.c,v 8.1 2015/04/20 11:14:20 hmb Exp $
 *****************************************************************************/

#include <stdlib.h>
#include <math.h>
#include <cfemdef.h>

#if defined(__uxp__)
#pragma global novrec
#pragma global noalias
#endif

#define MAX(x, y) ( ((x)>(y)) ? (x) : (y))


void dvamax (int_t n, 
	     const double* x, int_t incx,
	     const double* y, int_t incy,
	           double* z, int_t incz)
{
  register int_t i;
  register double  absx, absy;

  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;
  z += (incz<0) ? (-n+1)*incz : 0;

  for (i = 0; i < n; i++) {
    absx      = fabs (x[i*incx]);
    absy      = fabs (y[i*incy]);
    z[i*incz] = MAX (absx, absy);
  }
}


void ivamax (int_t n,
	     const int_t* x, int_t incx,
	     const int_t* y, int_t incy,
	           int_t* z, int_t incz)
{
  register int_t i, absx, absy;

  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;
  z += (incz<0) ? (-n+1)*incz : 0;

  for (i = 0; i < n; i++) {
    absx      = abs (x[i*incx]);
    absy      = abs (y[i*incy]);
    z[i*incz] = MAX (absx, absy);
  }
}


void svamax (int_t n,
	     const float* x, int_t incx,
	     const float* y, int_t incy,
	           float* z, int_t incz)
{
  register int_t i;
  register float   absx, absy;

  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;
  z += (incz<0) ? (-n+1)*incz : 0;

  for (i = 0; i < n; i++) {
#if  defined(__uxp__) || defined(_SX)
    absx = (float) fabs (x[i*incx]);
    absy = (float) fabs (y[i*incy]);
#else
    absx = fabsf (x[i*incx]);
    absy = fabsf (y[i*incy]);
#endif
    z[i*incz] = MAX (absx, absy);
  }
}
