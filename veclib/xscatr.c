/*****************************************************************************
 * xscatr:  vector scatter:  z[y[i]] = x[i].
 *
 * NB:  It is assumed that this operation is vectorizable, i.e. that there
 * are no repeated indices in the indirection vector y --- y is a permutator.
 *
 * $Id: xscatr.c,v 8.1 2015/04/20 11:14:20 hmb Exp $
 *****************************************************************************/

#include <cfemdef.h>

#if defined(__uxp__)
#pragma global novrec
#pragma global noalias
#endif
  
void dscatr (int_t n, const double* x, const int_t* y, double* z)
{
  register int_t i;

  for (i = 0; i < n; i++) z[y[i]] = x[i];
}


void iscatr (int_t n, const int_t* x, const int_t* y, int_t* z)
{
  register int_t i;

  for (i = 0; i < n; i++) z[y[i]] = x[i];
}


void sscatr (int_t n, const float* x, const int_t* y, float* z)
{
  register int_t i;

  for (i = 0; i < n; i++) z[y[i]] = x[i];
}
