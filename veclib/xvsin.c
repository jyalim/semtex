/*****************************************************************************
 * xvsin:  y[i] = sin(x[i]).
 *
 * $Id: xvsin.c,v 8.1 2015/04/20 11:14:21 hmb Exp $
 *****************************************************************************/

#include <math.h>
#include <cfemdef.h>

#if defined(__uxp__)
#pragma global novrec
#pragma global noalias
#endif

void dvsin (int_t n, const double* x, int_t incx,
                           double* y, int_t incy)
{
  register int_t i;

  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;

  for (i = 0; i < n; i++) y[i*incy] = sin (x[i*incx]);
}


void svsin (int_t n, const float* x, int_t incx,
                           float* y, int_t incy)
{
  register int_t i;

  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;

#if defined(__uxp__) || defined(_SX)
  for (i = 0; i < n; i++) y[i*incy] = (float) sin  (x[i*incx]);
#else
  for (i = 0; i < n; i++) y[i*incy] =         sinf (x[i*incx]);
#endif
}
