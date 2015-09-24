/*****************************************************************************
 * xramp:  ramp function:  x[i] = alpha + i*beta.
 *
 * $Id: xramp.c,v 8.1 2015/04/20 11:14:20 hmb Exp $
 *****************************************************************************/

#include <cfemdef.h>

#if defined(__uxp__)
#pragma global novrec
#pragma global noalias
#endif

void dramp (int_t n, double alpha, double beta, double* x, int_t incx)
{
  register int_t i;

  x += (incx<0) ? (-n+1)*incx : 0;

  for (i = 0; i < n; i++) x[i*incx] = alpha + i * beta;
}


void iramp (int_t n, int_t alpha, int_t beta, int_t* x, int_t incx)
{
  register int_t i;

  x += (incx<0) ? (-n+1)*incx : 0;

  for (i = 0; i < n; i++) x[i*incx] = alpha + i * beta;
}


void sramp (int_t n, float alpha, float beta, float* x, int_t incx)
{
  register int_t i;
  
  x += (incx<0) ? (-n+1)*incx : 0;

  for (i = 0; i < n; i++) x[i*incx] = alpha + i * beta;
}
