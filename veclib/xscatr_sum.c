/*****************************************************************************
 * xscatr_sum:  vector scatter with summation:  z[y[i]] += x[i].
 *
 * NB:  It is assumed that this operation is vectorizable, i.e. that there
 * are no repeated indices in the indirection vector y.
 *
 * $Id: xscatr_sum.c,v 8.1 2015/04/20 11:14:20 hmb Exp $
 *****************************************************************************/

#include <cfemdef.h>

#if defined(__uxp__)
#pragma global novrec
#pragma global noalias
#endif
  
void dscatr_sum (int_t n, const double* x, const int_t* y, double* z)
{
  register int_t i;

#if defined(__uxp__)
#pragma loop novrec z
#elif defined(_SX)
#pragma vdir nodep
#endif

  for (i = 0; i < n; i++) z[y[i]] += x[i];
}


void iscatr_sum (int_t n, const int_t* x, const int_t* y, int_t* z)
{
  register int_t i;

#if defined(__uxp__)
#pragma loop novrec z
#elif defined(_SX)
#pragma vdir nodep
#endif

  for (i = 0; i < n; i++) z[y[i]] += x[i];
}


void sscatr_sum (int_t n, const float* x, const int_t* y, float* z)
{
  register int_t i;

#if defined(__uxp__)
#pragma loop novrec z
#elif defined(_SX)
#pragma vdir nodep
#endif

  for (i = 0; i < n; i++) z[y[i]] += x[i];
}
