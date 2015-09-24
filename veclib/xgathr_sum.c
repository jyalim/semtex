/*****************************************************************************
 *  xgathr_sum:  vector gather, with summation:  z[i] += x[y[i]].
 *
 * $Id: xgathr_sum.c,v 8.1 2015/04/20 11:14:20 hmb Exp $
 *****************************************************************************/

#include <cfemdef.h>

#if defined(__uxp__)
#pragma global novrec
#pragma global noalias
#endif

void dgathr_sum (int_t n, const double* x, const int_t* y, double* z)
{
  register int_t i;

  for (i = 0; i < n; i++) z[i] += x[y[i]];
}


void igathr_sum (int_t n, const int_t* x, const int_t* y, int_t* z)
{
  register int_t i;

  for (i = 0; i < n; i++) z[i] += x[y[i]];
}


void sgathr_sum (int_t n, const float* x, const int_t* y, float* z)
{
  register int_t i;

  for (i = 0; i < n; i++) z[i] += x[y[i]];
}
