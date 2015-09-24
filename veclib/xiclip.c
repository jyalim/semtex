/*****************************************************************************
 * xiclip: inverted clip to interval [alpha,beta]:
 *   if x[i] < (alpha+beta)/2 y[i] = MIN(x[i],alpha) else y[i]=MAX(x[i],beta)
 *
 * $Id: xiclip.c,v 8.1 2015/04/20 11:14:20 hmb Exp $
 *****************************************************************************/

#include <cfemdef.h>

#if defined(__uxp__)
#pragma global novrec
#pragma global noalias
#endif

#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX(a, b) ((a) > (b) ? (a) : (b))

void diclip (int_t n, const double alpha, const double beta,
	     const double* x, int_t incx,
	           double* y, int_t incy)
{
  register int_t i;
  register double  xtmp;
  const double     mval = 0.5*(alpha + beta);

  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;

  for (i = 0; i < n; i++) {
    xtmp = x[i*incx];
    y[i*incy] = (xtmp < mval) ? MIN(xtmp, alpha) : MAX(xtmp, beta);
  }
}


void iiclip (int_t n, const int_t alpha, const int_t beta,
	     const int_t* x, int_t incx,
	           int_t* y, int_t incy)
{
  register int_t i, xtmp;
  const int_t    mval = (alpha + beta)/2;

  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;

  for (i = 0; i < n; i++) {
    xtmp = x[i*incx];
    y[i*incy] = (xtmp <= mval) ? MIN(xtmp, alpha) : MAX(xtmp, beta);
  }
}


void siclip (int_t n, const float alpha, const float beta,
	     const float* x, int_t incx,
	           float* y, int_t incy)
{
  register int_t i;
  register float   xtmp;
  const float      mval=0.5*(alpha + beta);

  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;

  for (i = 0; i < n; i++) {
    xtmp = x[i*incx];
    y[i*incy] = (xtmp < mval) ? MIN(xtmp, alpha) : MAX(xtmp, beta);
  }
}

#undef MIN
#undef MAX
