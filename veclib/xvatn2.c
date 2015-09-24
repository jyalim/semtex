/*****************************************************************************
 * xvatn2:  z[i] = atan2(x[i], y[i]).
 *
 * $Id: xvatn2.c,v 8.1 2015/04/20 11:14:20 hmb Exp $
 *****************************************************************************/

#include <math.h>
#include <cfemdef.h>

#if defined(__uxp__)
#pragma global novrec
#pragma global noalias
#endif

void dvatn2 (int_t n, 
	     const double* x, int_t incx,
	     const double* y, int_t incy,
	           double* z, int_t incz)
{
  register int_t i;

  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;
  z += (incz<0) ? (-n+1)*incz : 0;

  for (i = 0; i < n; i++) z[i*incz] = atan2 (x[i*incx], y[i*incy]);
}


void svatn2 (int_t n,
	     const float* x, int_t incx,
	     const float* y, int_t incy,
	           float* z, int_t incz)
{
  register int_t i;

  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;
  z += (incz<0) ? (-n+1)*incz : 0;

#if  defined(__uxp__) || defined(_SX)
  for (i = 0; i < n; i++) z[i*incz] = (float) atan2  (x[i*incx], y[i*incy]);
#else
  for (i = 0; i < n; i++) z[i*incz] =         atan2f (x[i*incx], y[i*incy]);
#endif
}
