/*****************************************************************************
 * Two- and three-component vector magnitudes.
 *
 * xvhypot:  z[i] = sqrt(SQR(x[i]) + SQR(y[i])).
 *
 * xvmag:    z[i] = sqrt(SQR(w[i]) + SQR(x[i]) + SQR(y[i])).
 *
 * $Id: xvhypot.c,v 8.1 2015/04/20 11:14:20 hmb Exp $
 *****************************************************************************/

#include <math.h>
#include <cfemdef.h>

#define SQR(a) ((a) * (a))

#if defined(__uxp__)
#pragma global novrec
#pragma global noalias
#endif

void dvhypot (int_t n,
	      const double* x, int_t incx,
	      const double* y, int_t incy,
	            double* z, int_t incz)
{
  register int_t i;

  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;
  z += (incz<0) ? (-n+1)*incz : 0;

  for (i = 0; i < n; i++) z[i*incz] = hypot (x[i*incx], y[i*incy]);
}


void svhypot (int_t n,
	      const float* x, int_t incx,
	      const float* y, int_t incy,
	            float* z, int_t incz)
{
  register int_t i;

  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;
  z += (incz<0) ? (-n+1)*incz : 0;

#if  defined (__uxp__) || defined(_SX)
  for (i = 0; i < n; i++) z[i*incz] = (float) hypot  (x[i*incx], y[i*incy]);
#else
  for (i = 0; i < n; i++) z[i*incz] =         hypotf (x[i*incx], y[i*incy]);
#endif
}


void dvmag (int_t n,
	    const double* w, int_t incw,
	    const double* x, int_t incx,
	    const double* y, int_t incy,
	          double* z, int_t incz)
{
  register int_t i;

  w += (incw<0) ? (-n+1)*incw : 0;
  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;
  z += (incz<0) ? (-n+1)*incz : 0;

  for (i = 0; i < n; i++)
    z[i*incz] = sqrt (SQR(w[i*incw]) + SQR(x[i*incx]) + SQR(y[i*incy]));
}


void svmag (int_t n,
	    const float* w, int_t incw,
	    const float* x, int_t incx,
	    const float* y, int_t incy,
	          float* z, int_t incz)
{
  register int_t i;

  w += (incw<0) ? (-n+1)*incw : 0;
  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;
  z += (incz<0) ? (-n+1)*incz : 0;

#if  defined (__uxp__) || defined(_SX)
  for (i = 0; i < n; i++)
    z[i*incz] = (float) sqrt  (SQR(w[i*incw])+SQR(x[i*incx])+SQR(y[i*incy]));
#else
  for (i = 0; i < n; i++)
    z[i*incz] =         sqrtf (SQR(w[i*incw])+SQR(x[i*incx])+SQR(y[i*incy]));
#endif
}

#undef SQR
