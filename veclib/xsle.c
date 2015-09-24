/*****************************************************************************
 * xsle:  y[i] = alpha <= x[i].
 *
 * $Id: xsle.c,v 8.1 2015/04/20 11:14:20 hmb Exp $
 *****************************************************************************/

#include <cfemdef.h>

#if defined(__uxp__)
#pragma global novrec
#pragma global noalias
#endif

void dsle (int_t n, double alpha,
	   const double* x, int_t incx,
	         double* y, int_t incy)
{
  register int_t i;

  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;

  for (i = 0; i < n; i++) y[i*incy] = (alpha <= x[i*incx]) ? 1 : 0;
}


void isle (int_t n, int_t alpha,
	   const int_t* x, int_t incx,
	         int_t* y, int_t incy)
{
  register int_t i;

  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;

  for (i = 0; i < n; i++) y[i*incy] = (alpha <= x[i*incx]) ? 1 : 0;
}


void ssle (int_t n, float alpha,
	   const float* x, int_t incx,
	         float* y, int_t incy)
{
  register int_t i;

  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;

  for (i = 0; i < n; i++) y[i*incy] = (alpha <= x[i*incx]) ? 1 : 0;
}
