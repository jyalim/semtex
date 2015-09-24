/*****************************************************************************
 * xvsin:  y[i] = tanh(x[i]).
 *
 * $Id: xvtanh.c,v 8.1 2015/04/20 11:14:21 hmb Exp $
 *****************************************************************************/

#include <math.h>
#include <cfemdef.h>

#if defined(__uxp__)
#pragma global novrec
#pragma global noalias
#endif

void dvtanh (int_t n, const double* x, int_t incx,
	                    double* y, int_t incy)
{
  register int_t i;

  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;

  for (i = 0; i < n; i++) y[i*incy] = tanh (x[i*incx]);
}


void svtanh (int_t n, const float* x, int_t incx,
                            float* y, int_t incy)
{
  register int_t i;

  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;

#if defined(__uxp__) || defined(_SX)
  for (i = 0; i < n; i++) y[i*incy] = (float) tanh  (x[i*incx]);
#else
  for (i = 0; i < n; i++) y[i*incy] =         tanhf (x[i*incx]);
#endif
}
