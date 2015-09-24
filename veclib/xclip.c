/*****************************************************************************
 * xclip: clip to interval [alpha,beta]: y[i] = MIN(MAX(x[i],alpha),beta).
 *
 * xclipup: clip on lower limit alpha:   y[i] = MAX(x[i],alpha).
 * xclipdn: clip on upper limit alpha:   y[i] = MIN(x[i],alpha).
 *
 * $Id: xclip.c,v 8.1 2015/04/20 11:14:20 hmb Exp $
 *****************************************************************************/

#include <cfemdef.h>

#if defined(__uxp__)
#pragma global novrec
#pragma global noalias
#endif

#define MIN(a, b)  ((a) < (b) ?     (a) : (b))
#define MAX(a, b)  ((a) > (b) ?     (a) : (b))

void dclip (int_t n,  const double alpha, const double beta,
	    const double* x, int_t incx,
	          double* y, int_t incy)
{
  register int_t i;

  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;

  for (i = 0; i < n; i++) y[i*incy] = MIN(MAX(x[i*incx],alpha),beta);
}


void iclip (int_t n,  const int_t alpha, const int_t beta,
	    const int_t* x, int_t incx,
	          int_t* y, int_t incy)
{
  register int_t i;

  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;

  for (i = 0; i < n; i++) y[i*incy] = MIN(MAX(x[i*incx],alpha),beta);
}


void sclip (int_t n,  const float alpha, const float beta,
	    const float* x, int_t incx,
	          float* y, int_t incy)
{
  register int_t i;

  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;

  for (i = 0; i < n; i++) y[i*incy] = MIN(MAX(x[i*incx],alpha),beta);
}


void dclipup (int_t n,  const double alpha,
	      const double* x, int_t incx,
	            double* y, int_t incy)
{
  register int_t i;

  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;

  for (i = 0; i < n; i++) y[i*incy] = MAX(x[i*incx], alpha);
}


void iclipup (int_t n,  const int_t alpha,
	      const int_t* x, int_t incx,
	            int_t* y, int_t incy)
{
  register int_t i;

  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;

  for (i = 0; i < n; i++) y[i*incy] = MAX(x[i*incx], alpha);
}


void sclipup (int_t n,  const float alpha,
	      const float* x, int_t incx,
	            float* y, int_t incy)
{
  register int_t i;

  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;

  for (i = 0; i < n; i++) y[i*incy] = MAX(x[i*incx], alpha);
}


void dclipdn (int_t n,  const double alpha,
	      const double* x, int_t incx,
	            double* y, int_t incy)
{
  register int_t i;

  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;

  for (i = 0; i < n; i++) y[i*incy] = MIN(x[i*incx], alpha);
}


void iclipdn (int_t n,  const int_t alpha,
	      const int_t* x, int_t incx,
	            int_t* y, int_t incy)
{
  register int_t i;

  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;

  for (i = 0; i < n; i++) y[i*incy] = MIN(x[i*incx], alpha);
}


void sclipdn (int_t n,  const float alpha,
	      const float* x, int_t incx,
	            float* y, int_t incy)
{
  register int_t i;

  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;

  for (i = 0; i < n; i++) y[i*incy] = MIN(x[i*incx], alpha);
}
