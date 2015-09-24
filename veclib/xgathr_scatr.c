/*****************************************************************************
 * xgathr_scatr:  vector gather-scatter:  z[y[i]] = w[x[i]].
 *
 * NB: it is assumed that this operation is vectorizable, i.e. that there
 * are no repeated indices in the indirection vector y[i].
 *
 * $Id: xgathr_scatr.c,v 8.1 2015/04/20 11:14:20 hmb Exp $
 *****************************************************************************/

#include <cfemdef.h>

#if defined(__uxp__)
#pragma global novrec
#pragma global noalias
#endif
  
void dgathr_scatr (int_t n,
		   const double* w,
		   const int_t* x, const int_t* y,
		   double* z)
{
  register int_t i;

  for (i = 0; i < n; i++) z[y[i]] = w[x[i]];
}


void igathr_scatr (int_t n,
		   const int_t* w,
		   const int_t* x, const int_t* y,
		   int_t* z)
{
  register int_t i;

  for (i = 0; i < n; i++) z[y[i]] = w[x[i]];
}


void sgathr_scatr (int_t n,
		   const float* w,
		   const int_t* x, const int_t* y,
		   float* z)
{
  register int_t i;

  for (i = 0; i < n; i++) z[y[i]] = w[x[i]];
}
