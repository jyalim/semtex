/*****************************************************************************
 * xvpow:  y[i] = pow(x[i], alpha).
 *
 * $Id: xspow.c,v 8.1 2015/04/20 11:14:20 hmb Exp $
 *****************************************************************************/

#include <math.h>
#include <cfemdef.h>

#if defined(__uxp__)
#pragma global novrec
#pragma global noalias
#endif

void dspow (const int_t n, const double alpha,
	    const double* x, int_t incx,
	          double* y, int_t incy)
{
  register int_t i;

  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;
  
  for (i = 0; i < n; i++) y[i*incy] = pow (x[i*incx], alpha);
}


void sspow (const int_t n, const float alpha,
	    const float* x, int_t incx, 
	          float* y, int_t incy)
{
  register int_t i;

  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;
  
  for (i = 0; i < n; i++)
#if  defined(__uxp__) || defined(_SX)
    y[i*incy] = (float) pow  (x[i*incx], alpha);
#else
    y[i*incy] =         powf (x[i*incx], alpha);
#endif
}
