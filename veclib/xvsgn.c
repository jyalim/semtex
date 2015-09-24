/*****************************************************************************
 * xvsgn:  y[i] = sign(x[i])  --  sign = -1 if x<0, else +1.
 *
 * $Id: xvsgn.c,v 8.1 2015/04/20 11:14:21 hmb Exp $
 *****************************************************************************/

#include <cfemdef.h>

#if defined(__uxp__)
#pragma global novrec
#pragma global noalias
#endif

void dvsgn (int_t n, const double* x, int_t incx,
                             double* y, int_t incy)
{
  register int_t i;

  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;

  for (i = 0; i < n; i++) y[i*incy] = (x[i*incx] < 0.0) ? -1.0 : 1.0;
}


void ivsgn (int_t n, const int_t* x, int_t incx,
                             int_t* y, int_t incy)
{
  register int_t i;

  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;

  for (i = 0; i < n; i++) y[i*incy] = (x[i*incx] < 0) ? -1 : 1;
}


void svsgn (int_t n, const float* x, int_t incx,
                             float* y, int_t incy)
{
  register int_t i;
  
  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;

  for (i = 0; i < n; i++) y[i*incy] = (x[i*incx] < 0.0F) ? -1.0F : 1.0F;
}
