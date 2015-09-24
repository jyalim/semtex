/*****************************************************************************
 * xvrecp:  y[i] = 1.0 / x[i].
 *
 * $Id: xvrecp.c,v 8.1 2015/04/20 11:14:21 hmb Exp $
 *****************************************************************************/

#include <cfemdef.h>

#if defined(__uxp__)
#pragma global novrec
#pragma global noalias
#endif

void dvrecp (int_t n, const double* x, int_t incx,
                              double* y, int_t incy)
{
  register int_t i;

  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;

  for (i = 0; i < n; i++) y[i*incy] = 1.0 / x[i*incx];
}


void svrecp (int_t n, const float* x, int_t incx,
                              float* y, int_t incy)
{
  register int_t i;

  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;

  for (i = 0; i < n; i++) y[i*incy] = 1.0F / x[i*incx];
}
