/*****************************************************************************
 * xvpow:  z[i] = pow(x[i], y[i]).
 *
 * $Id: xvpow.c,v 8.1 2015/04/20 11:14:21 hmb Exp $
 *****************************************************************************/

#include <math.h>
#include <cfemdef.h>

#if defined(__uxp__)
#pragma global novrec
#pragma global noalias
#endif

void dvpow (int_t n, 
	    const double* x, int_t incx,
	    const double* y, int_t incy,
	          double* z, int_t incz)
{
  register int_t i;

  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;
  z += (incz<0) ? (-n+1)*incz : 0;

  for (i = 0; i < n; i++) z[i*incz] = pow (x[i*incx], y[i*incy]);
}


void svpow (int_t n,
	    const float* x, int_t incx,
	    const float* y, int_t incy,
	          float* z, int_t incz)
{
  register int_t i;

  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;
  z += (incz<0) ? (-n+1)*incz : 0;

#if  defined(__uxp__) || defined(_SX)
  for (i = 0; i < n; i++) z[i*incz] = (float) pow  (x[i*incx], y[i*incy]);
#else
  for (i = 0; i < n; i++) z[i*incz] =         powf (x[i*incx], y[i*incy]);
#endif
}
