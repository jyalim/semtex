/*****************************************************************************
 * xznan:  x[i] = isnan(x[i]) ? 0.0 : x[i]
 *
 * $Id $
 *****************************************************************************/

#include <math.h>
#include <cfemdef.h>


#if defined(__uxp__)
#pragma global novrec
#pragma global noalias
#endif

void dznan (int_t n, double* x, int_t incx)
{
  register int_t i;

  x += (incx<0) ? (-n+1)*incx : 0;

#if defined(_SX)
#pragma vdir nodep
#endif

  for (i = 0; i < n; i++) x[i*incx] = isnan (x[i*incx]) ? 0.0 : x[i*incx];
}


void sznan (int_t n, float* x, int_t incx)
{
  register int_t i;

  x += (incx<0) ? (-n+1)*incx : 0;

#if defined(_SX)
#pragma vdir nodep
#endif

  for (i = 0; i < n; i++) x[i*incx] = isnan (x[i*incx]) ? 0.0 : x[i*incx];
}
